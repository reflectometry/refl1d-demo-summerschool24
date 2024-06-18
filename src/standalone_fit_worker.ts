// import { loadPyodideAndPackages } from './pyodide_worker.mjs';
import { expose } from 'comlink';
import { loadPyodide, version } from 'pyodide';
import type { PyodideInterface } from 'pyodide';
import { Signal } from './standalone_signal';
import type { PyProxy } from 'pyodide/ffi';
const DEBUG = true;


declare const REFL1D_WHEEL_FILE: string;
declare const BUMPS_WHEEL_FILE: string;
declare const MOLGROUPS_WHEEL_FILE: string;

type APIResult = PyProxy | number | string | boolean | null | undefined;


async function loadPyodideBase() {
    return await loadPyodide({
        indexURL: `https://cdn.jsdelivr.net/pyodide/v${version}/full/`
    });
}

async function loadBuiltins(pyodide: PyodideInterface) {
    await pyodide.loadPackage(["numpy", "scipy", "pytz", "h5py", "micropip"]);
}

async function doPipInstalls(pyodide: PyodideInterface) {
    await pyodide.runPythonAsync(`
    import micropip
    await micropip.install([
        "matplotlib",
        "periodictable",
        "blinker",
        "dill",
    ])
    `);
}

async function installLocalWheel(pyodide: PyodideInterface, wheel_file: string) {
    await pyodide.runPythonAsync(`
    import micropip
    await micropip.install("${wheel_file}", keep_going=True, deps=False)
    `);
}

async function createAPI(pyodide: PyodideInterface) {
    let api = await pyodide.runPythonAsync(`
    import dill
    from bumps.webview.server import api
    from refl1d.webview.server import api as refl1d_api
    api.state.parallel = 0
    api.state.problem.serializer = "dataclass"
    import refl1d

    # setup backend:
    refl1d.use('c_ext')

    # patch Thread to run in the main thread
    api.FitThread.start = api.FitThread.run
    api.FitThread.join = lambda self, timeout: None

    wrapped_api = {}

    def expose(method, method_name):
        def wrapper(args):
            pyargs = args.to_py() if args is not None else []
            return method(*pyargs)
        return wrapper

    for method_name, method in api.REGISTRY.items():
        if method_name in ["start_fit_thread"]:
            wrapped_api[method_name] = expose(method, method_name)

    def set_autosave_session_interval(interval: int):
        print("new interval", interval, type(interval))
        api.state.shared.autosave_session_interval = interval

    def set_problem(dilled_problem):
        problem = dill.loads(dilled_problem)
        api.state.problem.fitProblem = problem

    def get_nllf():
        return api.state.problem.fitProblem.nllf()

    wrapped_api["set_autosave_session_interval"] = expose(set_autosave_session_interval, "set_autosave_session_interval")
    wrapped_api["set_problem"] = expose(set_problem, "set_problem")
    wrapped_api["get_nllf"] = expose(get_nllf, "get_nllf")

    def fit_progress_handler(event):
        api.emit("evt_fit_progress", dill.dumps(event))

    def fit_complete_handler(event):
        api.emit("evt_fit_complete", dill.dumps(event))

    api.EVT_FIT_PROGRESS.connect(fit_progress_handler, weak=True)
    api.EVT_FIT_COMPLETE.connect(fit_complete_handler, weak=True)

    wrapped_api
`);
    return api;
}

async function loadPyodideAndPackages() { // loads pyodide
    const pyodide = await loadPyodideBase(); // run the function and wait for the result (base library)
    await loadBuiltins(pyodide); // waits until these python packpages are loaded to continue
    await doPipInstalls(pyodide);
    const preloaded_files_result = await fetch("../preloaded_files.json");
    let preloaded_files: { filename: string, path: string, source: string }[] = [];
    if (preloaded_files_result.ok) {
        preloaded_files = await preloaded_files_result.json() as { filename: string, path: string, source: string }[];
        for (let preload_file of preloaded_files) {
            if (preload_file.filename.endsWith(".whl")) {
                await installLocalWheel(pyodide, preload_file.source);
            }
        }
    }
    const api = await createAPI(pyodide);
    return api;
}

// export { loadPyodideAndPackages };

// let pyodideReadyPromise = loadPyodideAndPackages(); // run the functions stored in lines 4

type EventCallback = (message?: any) => any;

export class Server {
    handlers: { [signal: string]: EventCallback[] }
    pyodide: PyodideInterface;
    initialized: Promise<PyProxy>;
    nativefs: any;

    constructor() {
        this.handlers = {};
        this.nativefs = null;
        this.initialized = this.init();
    }

    async init() {
        const pyodide = await loadPyodideBase(); // run the function and wait for the result (base library)
        await loadBuiltins(pyodide); // waits until these python packpages are loaded to continue
        await doPipInstalls(pyodide);
        const preloaded_files_result = await fetch("../preloaded_files.json");
        let preloaded_files: { filename: string, path: string, source: string }[] = [];
        if (preloaded_files_result.ok) {
            preloaded_files = await preloaded_files_result.json() as { filename: string, path: string, source: string }[];
            for (let preload_file of preloaded_files) {
                if (preload_file.filename.endsWith(".whl")) {
                    await installLocalWheel(pyodide, preload_file.source);
                }
            }
        }
        const api = await createAPI(pyodide);
        this.pyodide = pyodide;
        const defineEmit = await pyodide.runPythonAsync(`
            def defineEmit(server):
                api.emit = server.asyncEmit;
            
            defineEmit
         `);
        await defineEmit(this);
        return api;
    }

    async set_signal(signal_in: Signal) {
        const api = await this.initialized;
        const { name, buffer } = signal_in;
        const signal = new Signal(name, buffer);
        console.log("setting abort signal in worker", signal);
        const defineFitEvent = await this.pyodide.runPythonAsync(`
            def defineFitEvent(event):
                api.state.${name} = event.to_py();
            
            defineFitEvent
        `);
        await defineFitEvent(signal);
    }

    async addHandler(signal: string, handler: EventCallback) {
        const signal_handlers = this.handlers[signal] ?? [];
        signal_handlers.push(handler);
        if (DEBUG) {
            console.log(`adding handler: ${signal}`);
        }
        if (signal === 'connect') {
            await this.initialized;
            await handler();
        }
        this.handlers[signal] = signal_handlers;
    }

    async removeHandler(signal: string, handler: EventCallback) {
        let signal_handlers = this.handlers[signal] ?? [];
        signal_handlers = signal_handlers.filter((h) => {
            if (h === handler) {
                console.log('matching worker handler found, removing: ', handler);
                return false;
            }
            return true;
        })
        this.handlers[signal] = signal_handlers;
    }

    async mount(dirHandle: FileSystemDirectoryHandle) {
        // const dirHandle = await self.showDirectoryPicker();
        console.log({dirHandle});   
        const nativefs = await this.pyodide.mountNativeFS("/home/pyodide/user_mount", dirHandle);
        this.nativefs = nativefs;
    }

    async syncFS() {
        let r = await this.nativefs?.syncfs?.();
    }

    async asyncEmit(signal: string, ...args: APIResult[]) {
        // this is for emit() calls from the python server
        const js_args = args.map((arg) => {
            return arg?.toJs?.({dict_converter: Object.fromEntries, create_pyproxies: false}) ?? arg;
        });
        const handlers = this.handlers[signal] ?? [];
        for (let handler of handlers) {
            handler(...js_args);
        }
    }

    async onAsyncEmit(signal: string, ...args: any[]) {
        // this is for responding to incoming emit() calls from the client
        // (which might be another server)
        const api = await this.initialized;
        const callback = (args[args.length - 1] instanceof Function) ? args.pop() : null;
        const result: PyProxy = await api.get(signal)(args);
        const jsResult = result?.toJs?.({dict_converter: Object.fromEntries, create_pyproxies: false}) ?? result;
        if (callback !== null) {
            await callback(jsResult);
        }
        result?.destroy?.();
        return jsResult;
    }

}

expose(Server);

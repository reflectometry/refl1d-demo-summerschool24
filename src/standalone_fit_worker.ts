// import { loadPyodideAndPackages } from './pyodide_worker.mjs';
import { expose } from 'comlink';
import { loadPyodide, version } from 'pyodide';
import type { PyodideInterface } from 'pyodide';
import { Signal } from './standalone_signal';
const DEBUG = true;

var pyodide: PyodideInterface;

declare const REFL1D_WHEEL_FILE: string;
declare const BUMPS_WHEEL_FILE: string;
declare const MOLGROUPS_WHEEL_FILE: string;

async function loadPyodideAndPackages() { // loads pyodide
    pyodide = await loadPyodide({
        indexURL: `https://cdn.jsdelivr.net/pyodide/v${version}/full/`
    }); // run the function and wait for the result (base library)

    await pyodide.loadPackage(["numpy", "scipy", "pytz", "h5py", "micropip"]); // waits until these python packpages are loaded to continue
  
    //import reductus library with micropip
    let api = await pyodide.runPythonAsync(`
    import micropip
    await micropip.install([
        "matplotlib",
        "plotly",
        "mpld3",
        "periodictable",
        "blinker",
        "dill",
    ])
    await micropip.install("../wheels/${BUMPS_WHEEL_FILE}")
    await micropip.install("../wheels/${REFL1D_WHEEL_FILE}", keep_going=True, deps=False)
    await micropip.install("../wheels/${MOLGROUPS_WHEEL_FILE}", keep_going=True, deps=False)

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

    def set_problem(dilled_problem):
        problem = dill.loads(dilled_problem)
        api.state.problem.fitProblem = problem

    wrapped_api["set_problem"] = expose(set_problem, "set_problem")

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

// export { loadPyodideAndPackages };

let pyodideReadyPromise = loadPyodideAndPackages(); // run the functions stored in lines 4

type EventCallback = (message?: any) => any;

export class Server {
    handlers: { [signal: string]: EventCallback[] }
    nativefs: any;

    constructor() {
        this.handlers = {};
        this.nativefs = null;
        this.init();
    }

    async init() {
        const api = await pyodideReadyPromise;
        const defineEmit = await pyodide.runPythonAsync(`
            def defineEmit(server):
                api.emit = server.asyncEmit;
            
            defineEmit
         `);
        await defineEmit(this);
    }

    async set_signal(signal_in: Signal) {
        const api = await pyodideReadyPromise;
        const { name, buffer } = signal_in;
        const signal = new Signal(name, buffer);
        console.log("setting abort signal in worker", signal);
        const defineFitEvent = await pyodide.runPythonAsync(`
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
            await pyodideReadyPromise;
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
        const nativefs = await pyodide.mountNativeFS("/home/pyodide/user_mount", dirHandle);
        this.nativefs = nativefs;
    }

    async syncFS() {
        let r = await this.nativefs?.syncfs?.();
    }

    async asyncEmit(signal: string, ...args: unknown[]) {
        // this is for emit() calls from the python server
        const js_args = args.map((arg) => {
            return arg?.toJs?.({dict_converter: Object.fromEntries}) ?? arg;
        });
        const handlers = this.handlers[signal] ?? [];
        for (let handler of handlers) {
            handler(...js_args);
        }
    }

    async onAsyncEmit(signal: string, ...args: any[]) {
        // this is for responding to incoming emit() calls from the client
        // (which might be another server)
        const api = await pyodideReadyPromise;
        const callback = (args[args.length - 1] instanceof Function) ? args.pop() : null;
        const result = await api.get(signal)(args);
        const jsResult = result?.toJs?.({dict_converter: Object.fromEntries}) ?? result;
        if (callback !== null) {
            await callback(jsResult);
        }
        return jsResult;
    }

}

expose(Server);

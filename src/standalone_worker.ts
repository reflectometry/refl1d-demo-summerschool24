// import { loadPyodideAndPackages } from './pyodide_worker.mjs';
import { expose, wrap, proxy, Remote } from 'comlink';
import { loadPyodide, version } from 'pyodide';
import type { PyodideInterface } from 'pyodide';
import type { PyProxy } from 'pyodide/ffi';
import { Signal } from './standalone_signal';
import type { Server as FitServer } from './standalone_fit_worker';

const DEBUG = true;

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
        "plotly",
        "mpld3",
        "periodictable",
        "blinker",
        "dill",
        "cloudpickle",
        "orsopy",
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
    import matplotlib
    matplotlib.use("agg")
    
    from typing import Any
    from bumps.webview.server import api
    import bumps.cli
    from refl1d.webview.server import api as refl1d_api
    from refl1d.bumps_interface import fitplugin
    api.state.parallel = 0
    api.state.problem.serializer = "dataclass"
    import molgroups
    print("api imported")

    import refl1d
    import asyncio
    import json
    import dill
    import cloudpickle

    # setup backend:
    bumps.cli.install_plugin(fitplugin)
    refl1d.use('c_ext')
    
    js_server_instance = object()
    wrapped_api = {}

    def expose(method, method_name):
        def wrapper(args):
            pyargs = args.to_py() if args is not None else []
            result = method(*pyargs)
            return result

        return wrapper

    for method_name, method in api.REGISTRY.items():
        print("wrapping:", method_name)
        wrapped_api[method_name] = expose(method, method_name)

    async def set_base_path(pathlist):
        # this is called from the client to set the base path
        # but unused in the pyodide environment
        pass

    async def worker_fit_progress_handler(serialized_event):
        event = dill.loads(serialized_event)
        await api._fit_progress_handler(event)
        if event.get("message") == "uncertainty_update":
            # call back into the JavaScript server to sync the filesystem
            await js_server_instance.syncFS()

    async def worker_fit_complete_handler(serialized_event):
        event = dill.loads(serialized_event)
        await api._fit_complete_handler(event)
        await js_server_instance.syncFS()

    wrapped_api["evt_fit_progress"] = expose(worker_fit_progress_handler, "evt_fit_progress")
    wrapped_api["evt_fit_complete"] = expose(worker_fit_complete_handler, "evt_fit_complete")
    wrapped_api["set_base_path"] = expose(set_base_path, "set_base_path")

    from dataclasses import dataclass

    @dataclass
    class WorkerFitThread:
        fitclass: Any
        fit_abort_event: Any
        problem: bumps.fitproblem.FitProblem
        mapper: Any
        options: dict
        parallel: int
        convergence_update: int
        uncertainty_update: int
        console_update: int = 500 # seconds
        fit_state: Any = None # for resuming fits
        convergence: Any = None # for resuming fits

        _alive = True

        def join(self, timeout=None):
            self._alive = False

        def is_alive(self):
            return self._alive

        def run(self):
            print("running dummy fit thread")
            asyncio.create_task(self._run())

        def start(self):
            print("started dummy fit thread")
            self.run()
        
        async def _run(self):
            dumped = cloudpickle.dumps(self.problem)
            await api.emit("set_fit_thread_autosave_session_interval", self.uncertainty_update)
            await api.emit("set_fit_thread_problem", dumped)
            await api.emit("start_fit_thread_fit", self.fitclass.id, self.options, False) # no resume for now
            await api.emit("add_notification", {
                "title": "Fit Started",
                "content": f"Fit started with problem: {api.state.problem.fitProblem.name}",
                "timeout": 2000,
            });

    api.FitThread = WorkerFitThread

    wrapped_api
    `);
    return api;
}

const fit_worker = new Worker(new URL("./standalone_fit_worker.ts", import.meta.url), {type: 'module'});
const FitServerClass = wrap<typeof FitServer>(fit_worker);
type EventCallback = (message?: any) => any;

export class Server {
    base_url: string;
    handlers: { [signal: string]: EventCallback[] }
    nativefs: any;
    pyodide: PyodideInterface;
    api: PyProxy;
    initialized: Promise<void>;
    fit_server_promise: Promise<Remote<FitServer>>;

    constructor(base_url: string) {
        this.base_url = base_url;
        this.handlers = {};
        this.nativefs = null;
        console.log(`starting server at ${base_url}`);
        this.fit_server_promise = new FitServerClass(base_url);
        this.initialized = this.init();
    }

    async init() {
        await this.asyncEmit("server_startup_status", {status: "loading python", percent: 0});
        const pyodide = await loadPyodideBase();
        this.pyodide = pyodide;
        
        await this.asyncEmit("server_startup_status", {status: "initializing builtin modules", percent: 25});
        await loadBuiltins(pyodide);
        await this.asyncEmit("server_startup_status", {status: "installing pip dependencies", percent: 50});
        await doPipInstalls(pyodide);
        const preloaded_files_result = await fetch("../preloaded_files.json");
        if (preloaded_files_result.ok) {
            const preloaded_files = await preloaded_files_result.json() as { filename: string, path: string, source: string }[];
            for (let preload_file of preloaded_files) {
                if (preload_file.path.startsWith("public/examples")) {
                    const target_path = preload_file.path.replace(/^public\/examples/, "/home/pyodide");
                    await pyodide.FS.mkdirTree(target_path);
                    await pyodide.FS.createLazyFile(target_path, preload_file.filename, preload_file.source, true, false);
                }
                if (preload_file.filename.endsWith(".whl")) {
                    console.log("installing wheel:", `${this.base_url}${preload_file.source}`);
                    await installLocalWheel(pyodide, `${this.base_url}${preload_file.source}`);
                }
                
            }
        }
        await this.asyncEmit("server_startup_status", {status: "loading bumps and refl1d", percent: 75});
        const api = await createAPI(pyodide);
        this.api = api;
        await this.asyncEmit("server_startup_status", {status: "api created", percent: 100});
        const fit_server = await this.fit_server_promise;
        const abort_fit_signal = new Signal("fit_abort_event");
        const fit_complete_signal = new Signal("fit_complete_event");
        await fit_server.set_signal(abort_fit_signal);
        await this.set_signal(abort_fit_signal);
        await fit_server.set_signal(fit_complete_signal);
        await this.set_signal(fit_complete_signal);
        const defineJSServer = await pyodide.runPythonAsync(`
            def defineJSServer(js_server):
                global js_server_instance
                js_server_instance = js_server
                api.emit = js_server.asyncEmit;
            
            defineJSServer
        `);
        await defineJSServer(this);
        this.addHandler('set_fit_thread_autosave_session_interval', async (interval: number) => {
            const result = await fit_server.onAsyncEmit('set_autosave_session_interval', interval);
        });
        this.addHandler('set_fit_thread_problem', async (problem: any) => {
            const result = await fit_server.onAsyncEmit('set_problem', problem);
            console.log("set_fit_thread_problem result:", result);
        });
        this.addHandler('start_fit_thread_fit', async (...args: any[]) => {
            const result = await fit_server.onAsyncEmit('start_fit_thread', ...args);
            console.log("start_fit_thread_fit result:", result);
        });
        const fit_progress_handler = async (event: any) => {
            await this.onAsyncEmit('evt_fit_progress', event);
        }
        fit_server.addHandler('evt_fit_progress', proxy(fit_progress_handler));
        const fit_complete_handler = async (event: any) => {
            await this.onAsyncEmit('evt_fit_complete', event);
        }
        fit_server.addHandler('evt_fit_complete', proxy(fit_complete_handler));
    }

    async set_signal(signal_in: Signal) {
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
            await this.fit_server_promise;
            await handler();
        }
        this.handlers[signal] = signal_handlers;
    }

    async removeHandler(signal: string, handler: EventCallback) {
        let signal_handlers = this.handlers[signal] ?? [];
        signal_handlers = signal_handlers.filter((h) => {
            if (h === handler) {
                console.log('matching worker handler found, removing: ', handler);
                handler?.destroy?.();
                return false;
            }
            return true;
        })
        this.handlers[signal] = signal_handlers;
    }

    async mount(dirHandle: FileSystemDirectoryHandle) {
        const { FS } = this.pyodide;
        const mountpoint = "/home/pyodide/local_mount";
        const analyze = FS.analyzePath(mountpoint);
        if (analyze.exists && analyze?.object?.mount?.mountpoint === mountpoint) {
            console.log("unmounting existing mountpoint");
            FS.unmount(mountpoint);
        }
        const nativefs = await this.pyodide.mountNativeFS(mountpoint, dirHandle);
        this.nativefs = nativefs;
        this.asyncEmit("base_path", mountpoint.split("/").map((part) => (part === '') ? '/' : part));
    }

    async syncFS() {
        let r = await this.nativefs?.syncfs?.();
    }

    async asyncEmit(signal: string, ...args: unknown[]) {
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
        // this is for emit() calls from the client
        await this.initialized; // api ready after this...
        const callback = (args[args.length - 1] instanceof Function) ? args.pop() : null;
        console.debug(`onAsyncEmit: ${signal}`, {args, callback});
        const api_function = this.api.get(signal);
        if (api_function === undefined) {
            throw new Error(`Signal ${signal} not found in API`);
        }
        else if (typeof api_function !== 'function') {
            throw new Error(`Signal ${signal} is not a function in API`);
        }
        console.debug("api function:", api_function, "args:", args);

        const result: PyProxy | undefined | number | string | null = await api_function(args);
        const jsResult = result?.toJs?.({dict_converter: Object.fromEntries, create_pyproxies: false}) ?? result;
        if (callback !== null) {
            await callback(jsResult);
        }
        result?.destroy?.();
        return jsResult;
    }

}

expose(Server);

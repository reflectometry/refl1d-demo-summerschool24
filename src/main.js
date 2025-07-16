import "bootstrap/dist/css/bootstrap.min.css";
import { computed, createApp } from "vue";
import { file_menu_items, fileBrowser, shared_state } from "bumps-webview-client/src/app_state";
import App from "bumps-webview-client/src/App.vue";
import { io } from "./asyncWorkerSocket";
import { dqIsFWHM } from "refl1d-webview-client/src/app_state";
import { panels } from "refl1d-webview-client/src/panels";
import "bumps-webview-client/src/style.css";

const name = "Refl1D";
const urlParams = new URLSearchParams(window.location.search);
const singlePanel = urlParams.get("single_panel");

const can_mount_local = ( 'showDirectoryPicker' in window );

const socket = io();

async function mountLocal() {
  const success = (await socket.mountLocal?.()) ?? false;
}

async function loadProbeFromFile() {
  if (fileBrowser.value) {
    const settings = {
      title: "Load Probe Data from File",
      callback: (pathlist, filename) => {
        socket.asyncEmit("load_probe_from_file", pathlist, filename, 0, dqIsFWHM.value);
      },
      show_name_input: true,
      name_input_label: "Filename",
      require_name: true,
      show_files: true,
      chosenfile_in: "",
      search_patterns: [".ort", ".orb", ".refl", ".dat", ".txt"],
    };
    fileBrowser.value.open(settings);
  }
}

createApp(App, { panels, socket, singlePanel, name }).mount("#app");

const modelNotLoaded = computed(() => shared_state.model_file == null);

const quit_item_index = file_menu_items.value.findIndex(menu_item => menu_item.text === "Quit");
if (quit_item_index !== -1) {
  file_menu_items.value.splice(quit_item_index, 1);
}

file_menu_items.value.push({
  text: "Load Data into Model",
  action: loadProbeFromFile,
  disabled: modelNotLoaded,
});
file_menu_items.value.push({
    text: "Mount Local Directory",
    action: mountLocal,
    disabled: !can_mount_local,
    icon: "folder-open",
    tooltip: "Mount a local directory to access files in the Refl1D webview",
});

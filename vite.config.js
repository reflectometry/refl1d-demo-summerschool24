import { defineConfig } from "vite";
import svgLoader from "vite-svg-loader";
import vue from "@vitejs/plugin-vue";

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [vue(), svgLoader()],
  base: "",
  define: {
    // By default, Vite doesn't include shims for NodeJS/
    // necessary for segment analytics lib to work
    global: {},
    APP_VERSION: JSON.stringify(process.env.npm_package_version),
  },
  worker: {
    format: 'es',
    rollupOptions: {
      external: ["node-fetch"],
    },
  },
});
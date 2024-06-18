import { fileURLToPath, URL } from 'node:url'
import { join } from 'node:path'

import { defineConfig } from 'vite'
import vue from '@vitejs/plugin-vue'

// https://vitejs.dev/config/
export default ({mode}) => {
  return defineConfig({
    plugins: [
      vue(),
    ],
    base: '',
    resolve: {
      alias: {
        '@': fileURLToPath(new URL('./src', import.meta.url)),
        './asyncSocket.ts': join(process.cwd(), 'src', 'asyncWorkerSocket.ts'),
        'socket.io-client': join(process.cwd(), 'src', 'asyncWorkerSocket.ts'),
      }
    },
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
    build: {
      rollupOptions: {
        output: {
          // Default
          // dir: (mode == 'standalone') ? 'dist_standalone': join('dist', process.env.npm_package_version),
          entryFileNames: (mode == 'production') ? 'assets/[name].js' : 'assets/[name].[hash].js',
          assetFileNames: (mode == 'production') ? 'assets/[name][extname]' : 'assets/[name].[hash][extname]',
          // chunkFileNames: "chunk-[name].js",
          // manualChunks: undefined,
        }
      }
    }
  })
}

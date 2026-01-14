const CACHE_NAME = 'eos-solver-v1';

// We cache the HTML, the Tailwind CSS (optional if you want offline styles), 
// and the critical Pyodide files.
const urlsToCache = [
  './',
  './index.html',
  './manifest.json',
  'https://cdn.tailwindcss.com',
  'https://cdn.jsdelivr.net/pyodide/v0.25.0/full/pyodide.js',
  'https://cdn.jsdelivr.net/pyodide/v0.25.0/full/pyodide.asm.js',
  'https://cdn.jsdelivr.net/pyodide/v0.25.0/full/pyodide.asm.wasm',
  'https://cdn.jsdelivr.net/pyodide/v0.25.0/full/repodata.json'
];

self.addEventListener('install', event => {
  event.waitUntil(
    caches.open(CACHE_NAME)
      .then(cache => cache.addAll(urlsToCache))
  );
});

self.addEventListener('fetch', event => {
  event.respondWith(
    caches.match(event.request)
      .then(response => {
        // Return cached file if found, otherwise fetch from network
        return response || fetch(event.request);
      })
  );
});

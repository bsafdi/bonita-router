/* Bonita Router service worker -- cache the shell so the app opens with no signal. */
const CACHE = 'bonita-20260912154136';
const ASSETS = ['./', './index.html', './manifest.webmanifest', './icon.svg'];
self.addEventListener('install', e => {
  self.skipWaiting();
  e.waitUntil(caches.open(CACHE).then(c => c.addAll(ASSETS)).catch(()=>{}));
});
self.addEventListener('activate', e => {
  e.waitUntil(caches.keys().then(ks =>
    Promise.all(ks.filter(k => k !== CACHE).map(k => caches.delete(k)))).then(()=>self.clients.claim()));
});
self.addEventListener('fetch', e => {
  const u = new URL(e.request.url);
  if (e.request.method !== 'GET') return;
  /* NOAA / Open-Meteo: network first, fall back to the last good response */
  if (u.origin !== location.origin) {
    e.respondWith(fetch(e.request).then(r => {
      const c = r.clone(); caches.open(CACHE).then(k => k.put(e.request, c));
      return r;
    }).catch(() => caches.match(e.request)));
    return;
  }
  /* the app shell: cache first, so it opens instantly and offline */
  e.respondWith(caches.match(e.request).then(r => r || fetch(e.request)));
});

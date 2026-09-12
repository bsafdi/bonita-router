/* Bonita Router service worker -- cache the shell so the app opens with no signal. */
const CACHE = 'bonita-20260912161958';
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
  /* The app shell: serve the cache immediately so it opens instantly and with
     no signal, but ALWAYS refetch in the background and store the result, so
     the next open is current.  Pure cache-first stranded phones on a stale
     build: nothing refetched the shell until the browser happened to re-check
     sw.js, which an installed PWA does lazily. */
  e.respondWith(caches.match(e.request).then(r => {
    const net = fetch(e.request).then(res => {
      if (res && res.ok) { const c = res.clone(); caches.open(CACHE).then(k => k.put(e.request, c)); }
      return res;
    }).catch(() => r);
    return r || net;
  }));
});

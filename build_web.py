#!/usr/bin/env python3
"""
Build the two shipping forms of the on-water app from one source.

  dist/index.html   complete standalone document for your own https host
                    (GitHub Pages, Netlify, anything).  Live NOAA + Open-Meteo
                    fetch, geolocation, installable PWA with offline cache.
  web/artifact.html content-only form for publishing as a Claude artifact:
                    same app, but that sandbox blocks outbound fetch, so it
                    runs on the NOAA predictions baked into model_data.json.

Both are single files -- model.js and model_data.json are inlined.
"""
import json, os, pathlib, datetime, subprocess, sys

HERE = pathlib.Path(__file__).parent

# Regenerate the browser app's data blob from the Python package, so the two
# implementations cannot drift apart.
subprocess.run([sys.executable, str(HERE / "export_model_data.py")], check=True)
body = (HERE / "web/app.body.html").read_text()
model = (HERE / "web/model.js").read_text()
data = json.loads((HERE / "web/model_data.json").read_text())

inline = ("<script>\n" + model + "\nwindow.BONITA_DATA = "
          + json.dumps(data, separators=(",", ":")) + ";\n</script>")
body = body.replace('<script src="MODEL_JS_PLACEHOLDER"></script>', inline)

(HERE / "web/artifact.html").write_text(body)

dist = HERE / "dist"
dist.mkdir(exist_ok=True)
doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">
<meta name="theme-color" content="#0b1220">
<meta name="description" content="Time-optimal routing for the BYC Big Windward-Leeward Race to Point Bonita, with live NOAA currents and wind.">
<meta name="apple-mobile-web-app-capable" content="yes">
<meta name="apple-mobile-web-app-status-bar-style" content="black-translucent">
<link rel="manifest" href="manifest.webmanifest">
<style>html,body{{margin:0;padding:0;background:#0b1220;color:#e8eef6;
font:14px/1.4 system-ui,sans-serif}}img{{max-width:100%}}[hidden]{{display:none!important}}</style>
{body}
<script>
if ('serviceWorker' in navigator && location.protocol === 'https:')
  window.addEventListener('load', () => navigator.serviceWorker.register('sw.js').catch(()=>{{}}));
</script>
</body>
</html>
"""
(dist / "index.html").write_text(doc)

(dist / "manifest.webmanifest").write_text(json.dumps({
    "name": "Bonita Router", "short_name": "Bonita",
    "description": "Time-optimal routing for the BYC Big Windward-Leeward to Point Bonita",
    "start_url": ".", "display": "standalone", "orientation": "portrait",
    "background_color": "#0b1220", "theme_color": "#0b1220",
    "icons": [{"src": "icon.svg", "sizes": "any", "type": "image/svg+xml", "purpose": "any"}],
}, indent=1))

(dist / "icon.svg").write_text(
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 192 192">'
    '<rect width="192" height="192" rx="30" fill="#0b1220"/>'
    '<path d="M96 26 L138 132 H54 Z" fill="#f2b233"/>'
    '<path d="M30 148 q33 14 66 0 t66 0" stroke="#4fa8e8" stroke-width="9" fill="none" '
    'stroke-linecap="round"/></svg>')

ver = datetime.datetime.utcnow().strftime("%Y%m%d%H%M%S")
(dist / "sw.js").write_text(f"""/* Bonita Router service worker -- cache the shell so the app opens with no signal. */
const CACHE = 'bonita-{ver}';
const ASSETS = ['./', './index.html', './manifest.webmanifest', './icon.svg'];
self.addEventListener('install', e => {{
  self.skipWaiting();
  e.waitUntil(caches.open(CACHE).then(c => c.addAll(ASSETS)).catch(()=>{{}}));
}});
self.addEventListener('activate', e => {{
  e.waitUntil(caches.keys().then(ks =>
    Promise.all(ks.filter(k => k !== CACHE).map(k => caches.delete(k)))).then(()=>self.clients.claim()));
}});
self.addEventListener('fetch', e => {{
  const u = new URL(e.request.url);
  if (e.request.method !== 'GET') return;
  /* NOAA / Open-Meteo: network first, fall back to the last good response */
  if (u.origin !== location.origin) {{
    e.respondWith(fetch(e.request).then(r => {{
      const c = r.clone(); caches.open(CACHE).then(k => k.put(e.request, c));
      return r;
    }}).catch(() => caches.match(e.request)));
    return;
  }}
  /* The app shell: serve the cache immediately so it opens instantly and with
     no signal, but ALWAYS refetch in the background and store the result, so
     the next open is current.  Pure cache-first stranded phones on a stale
     build: nothing refetched the shell until the browser happened to re-check
     sw.js, which an installed PWA does lazily. */
  e.respondWith(caches.match(e.request).then(r => {{
    const net = fetch(e.request).then(res => {{
      if (res && res.ok) {{ const c = res.clone(); caches.open(CACHE).then(k => k.put(e.request, c)); }}
      return res;
    }}).catch(() => r);
    return r || net;
  }}));
}});
""")

(dist / ".nojekyll").write_text("")
for f in sorted(dist.iterdir()):
    print(f"  {f.name:24s} {f.stat().st_size/1024:8.1f} KB")
print(f"  web/artifact.html        {(HERE/'web/artifact.html').stat().st_size/1024:8.1f} KB")

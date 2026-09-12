# Deploying the Bonita Router

Everything in this folder is static. There is no server, no build step, no key.
The page fetches NOAA and Open-Meteo directly from your phone — I verified both
allow cross-origin browser requests.

## GitHub Pages (recommended, ~5 minutes)

```bash
cd dist
git init && git add -A && git commit -m "Bonita router"
gh repo create bonita-router --public --source=. --push     # or create it in the web UI
gh api -X POST repos/:owner/bonita-router/pages -f source[branch]=main -f source[path]=/
```

Then open `https://<you>.github.io/bonita-router/` on your phone and use
**Share → Add to Home Screen**. That installs it as an app: full screen, its own
icon, and the service worker keeps the shell cached so it opens with no signal.

Netlify or Cloudflare Pages work identically — drag this folder onto their
dashboard. Any static https host is fine. It must be **https**, not `file://` —
browsers only grant geolocation on a secure origin.

## What it fetches, live

| Source | What | Refresh |
|---|---|---|
| NOAA CO-OPS `currents_predictions` | Harmonic current predictions at 8 stations in the race area | Once per day, cached |
| NOAA CO-OPS `wind` | 6-minute observed wind at San Francisco (37.806, −122.466 — at the Gate), Richmond, Alameda | Each refresh |
| Open-Meteo | Gridded wind forecast over 6 points across the Bay | Once, cached |

If a fetch fails the app falls back to `localStorage`, then to the NOAA
predictions baked into the page at build time. The status chip in the top bar
says which: **LIVE**, **CACHED**, or **BAKED**.

## SFBOFS (optional, the one thing a browser can't do)

NOAA's hydrodynamic current model is OPeNDAP/netCDF — not fetchable from a page.
Race morning, run this on your laptop and commit the result:

```bash
python prefetch_sfbofs.py --date 2026-09-12 --out dist/sfbofs.json
git -C dist commit -am "sfbofs $(date +%F)" && git -C dist push
```

The app picks up `sfbofs.json` if it's there. Skip it and the 8-station harmonic
model runs instead, which is what the whole analysis has used to-date.

## On the water

- **Race state** — switch Pre-start → To Bonita → To finish as you go.
- **Position** — *Start GPS* once, then leave it. *Simulate* flies the boat down
  the planned route so you can rehearse the display on shore.
- **Wind** — leave on *Live + forecast*. If you can see the breeze is different
  from what it says, switch to *Manual* and set what you actually have; the app
  divides your reading by the field's structure where you are, so the Gate
  acceleration and the island lees still apply on top.
- **Calibration** — the app compares your GPS speed against what it predicts for
  your course and position. If the ratio sits away from 1.00 for a few minutes,
  tap *Apply as wind correction*. That is usually the honest reading; the polar
  is the other candidate and there's a slider for it.
- It re-routes every 20 s from wherever you are, and on every leg change.

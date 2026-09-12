# Getting this onto GitHub Pages

You need the app at an **https** URL — browsers only grant geolocation on a
secure origin, and `file://` will not do.

```bash
unzip bonita_router.zip && cd bonita
git init && git add -A && git commit -m "Bonita router"
gh repo create bonita-router --private --source=. --push
```

Then in the repo: **Settings → Pages → Source: GitHub Actions**. The included
workflow rebuilds `dist/` from the Python tables and publishes on every push to
`main`, so the site can never drift from the model.

Open `https://<you>.github.io/bonita-router/` on your phone and use
**Share → Add to Home Screen**.

If you would rather not wait on Actions, `dist/` is already built — drag that
folder onto Netlify or Cloudflare Pages and you are done in a minute.

## Then, with Claude Code in the repo

`CLAUDE.md` carries the architecture, the traps and the known-weak inputs, so
Claude Code picks them up automatically. Useful first asks:

- "Put the signal boat's real position into `MARKS['start']` and re-run both courses."
- "Refetch tomorrow's NOAA currents and rebuild."  (`rm data/cache/currents_*` then `python run_race.py`)
- "Widen the Cityfront eddy in `DERIVED_ZONES` to a band of 5 seeds from Fort Point to Fort Mason and tell me what it changes."
- "Pull surface-bin current predictions for SFB1201/1202/1203 and see if the ebb asymmetry goes away."

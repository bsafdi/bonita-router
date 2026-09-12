"""
Geometry for the Central San Francisco Bay: coastline, water mask, projected grid.

The coastline here is HAND-DIGITIZED at roughly 100-300 m fidelity from chart
knowledge.  It is good enough for routing on a 200-300 m grid, but it is not a
survey.  If you want to do better, replace `LAND_POLYGONS` with polygons built
from NOAA ENC US5CA12M / OSM coastline; everything downstream only needs a
callable land mask.

Coordinates are (lon, lat) degrees.  The projection is a local equirectangular
approximation about LAT0 -- over a 25 km domain the distortion is < 0.1%.
"""

from __future__ import annotations

import numpy as np

# --------------------------------------------------------------------------
# Domain and projection
# --------------------------------------------------------------------------

LAT0 = 37.85          # projection reference latitude
LON0 = -122.42        # projection reference longitude

M_PER_DEG_LAT = 111_132.0
M_PER_DEG_LON = 111_320.0 * np.cos(np.deg2rad(LAT0))

# Bounding box of the model domain.  West edge is offshore of Point Bonita so
# the mark rounding has sea room; east edge is the Berkeley shoreline.
DOMAIN = dict(lon_min=-122.570, lon_max=-122.295, lat_min=37.770, lat_max=37.920)


def ll_to_xy(lon, lat):
    """Lon/lat degrees -> local metres (x east, y north)."""
    x = (np.asarray(lon) - LON0) * M_PER_DEG_LON
    y = (np.asarray(lat) - LAT0) * M_PER_DEG_LAT
    return x, y


def xy_to_ll(x, y):
    lon = np.asarray(x) / M_PER_DEG_LON + LON0
    lat = np.asarray(y) / M_PER_DEG_LAT + LAT0
    return lon, lat


# --------------------------------------------------------------------------
# Course marks  (edit these if the SIs put them elsewhere)
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Course marks, from the Sailing Instructions for 12 September 2026.
#
#   SI 13.0   "PB"        Point Bonita RG Fl (2+1) R 6s Bell
#                         37 48.70'N / 122 32.40'W
#             YRA 15      "EASOM", large lighted yellow column buoy off
#             "EASOM"     Yellow Bluff west of Harding Rock
#                         37 50.49'N / 122 28.08'W
#   SI 12.3   ALL MARKS ARE ROUNDED TO PORT.
#   SI 12.1   The course is signalled by numeral pennant at the warning
#             signal, so you do not know which one you are sailing until
#             11:00.  Both are modelled; pick with --course.
#   SI 5.1    Signal boat anchors on the Berkeley Circle about 1.5 NM west
#             of the Berkeley Marina breakwater, between YRA DOC and FOC.
#             VERIFY the actual position on the day -- this is a best
#             estimate from the SI text and Figure 1.
# --------------------------------------------------------------------------

# YRA Berkeley Circle marks that bracket the start (SI 5.1):
#   DOC  = "Berkeley Marina Channel Light 3", Fl G 4s 15ft 4M "3".  NOAA ENC
#          US5OAKFH light feature at 37.85767N 122.34943W; the YRA mark list
#          (github nkiesel/YRA-Marks) gives the same position.  VERIFIED.
#   FOC  = Olympic Circle buoy "F", 37.865N 122.37667W from the YRA mark list
#          (4-decimal precision, about +-100 m).  Not a charted aid.
YRA_DOC = (-122.34943, 37.85767)
YRA_FOC = (-122.37667, 37.86500)

MARKS = {
    # SI 5.1: signal boat ~1.5 NM W of the Berkeley Marina breakwater, between
    # DOC and FOC.  DOC itself is 1.5 NM from the breakwater end (ENC light at
    # 37.8672N 122.3213W), FOC is 2.6 NM, so the boat is near the DOC end of
    # that line.  This point is on the DOC-FOC line 500 m from DOC (1.66 NM
    # from the breakwater).  ESTIMATE, +-0.5 NM along the DOC-FOC line.
    "start":  (-122.3548, 37.8591),
    "finish": (-122.3548, 37.8591),      # SI 6.1, same line as the start
    "bonita": (-122.540000, 37.811667),  # SI 13.0 "PB"   37 48.70'N 122 32.40'W
    "easom":  (-122.468000, 37.841500),  # SI 13.0 YRA 15 37 50.49'N 122 28.08'W
}

# SI 12.0: Course 1 rounds "PB" (~20 NM), Course 2 rounds "EASOM" (~15 NM).
COURSES = {
    1: {"mark": "bonita", "name": 'YRA "PB" (Point Bonita)', "approx_nm": 20},
    2: {"mark": "easom",  "name": 'YRA 15 "EASOM" (Yellow Bluff)', "approx_nm": 15},
}

# SI 3.1 first warning 1100; the gun follows the warning by five minutes.
FIRST_WARNING_H = 11.0
DEFAULT_START_H = 11.0 + 5.0 / 60.0
TIME_LIMIT_H = 19.0                      # SI 6.2: finish after 1900 = DNF

MARK_RADIUS_M = 200.0   # how close counts as "rounded"


# --------------------------------------------------------------------------
# Land polygons (lon, lat), each a closed ring.  Inland filler vertices are
# deliberately crude -- only the waterline matters.
# --------------------------------------------------------------------------

MARIN = [
    # Outer (Pacific) coast, north to south
    (-122.610, 37.920), (-122.560, 37.905), (-122.545, 37.890),
    (-122.537, 37.872), (-122.535, 37.860), (-122.539, 37.845),
    (-122.539, 37.832),                                   # Rodeo Beach
    (-122.5294, 37.8158),                                 # Point Bonita
    (-122.522, 37.8175), (-122.514, 37.820),              # Bonita Cove
    (-122.5010, 37.8260),                                 # Point Diablo
    (-122.4930, 37.8272), (-122.4890, 37.8290),           # Kirby Cove
    (-122.4830, 37.8268),
    (-122.4785, 37.8255),                                 # Lime Point / north tower
    (-122.4770, 37.8300),
    (-122.4760, 37.8380),                                 # Yellow Bluff
    (-122.4775, 37.8430),
    (-122.4790, 37.8480),                                 # Sausalito waterfront
    (-122.4805, 37.8530),
    (-122.4830, 37.8600),
    (-122.4855, 37.8650),                                 # Richardson Bay entrance (N side)
    (-122.4810, 37.8690),
    (-122.4740, 37.8700),
    (-122.4690, 37.8715),                                 # Belvedere Point
    (-122.4645, 37.8705),
    (-122.4620, 37.8710),                                 # Peninsula Point
    (-122.4600, 37.8735),
    (-122.4560, 37.8740),                                 # Tiburon
    (-122.4500, 37.8760),
    (-122.4460, 37.8800),
    (-122.4530, 37.8880),                                 # Paradise Cay side
    (-122.4600, 37.8960),
    (-122.4700, 37.9100),
    # inland filler
    (-122.560, 37.940), (-122.640, 37.940), (-122.640, 37.900),
]

ANGEL_ISLAND = [
    (-122.4460, 37.8620),   # Point Stuart
    (-122.4430, 37.8570),
    (-122.4400, 37.8540),   # Point Knox
    (-122.4340, 37.8515),
    (-122.4270, 37.8505),
    (-122.4200, 37.8510),   # Point Blunt
    (-122.4165, 37.8555),
    (-122.4180, 37.8600),   # Point Simpton
    (-122.4200, 37.8655),
    (-122.4240, 37.8700),   # Point Campbell
    (-122.4310, 37.8710),
    (-122.4380, 37.8680),
    (-122.4400, 37.8690),   # Point Ione
    (-122.4440, 37.8660),
]

ALCATRAZ = [
    (-122.4265, 37.8280), (-122.4225, 37.8285), (-122.4195, 37.8270),
    (-122.4210, 37.8252), (-122.4255, 37.8252),
]

TREASURE_ISLAND = [
    (-122.3770, 37.8180), (-122.3745, 37.8300), (-122.3700, 37.8325),
    (-122.3640, 37.8300), (-122.3628, 37.8195), (-122.3745, 37.8175),
]

YERBA_BUENA = [
    (-122.3745, 37.8175), (-122.3628, 37.8195), (-122.3600, 37.8130),
    (-122.3620, 37.8055), (-122.3690, 37.8035), (-122.3730, 37.8090),
]

SAN_FRANCISCO = [
    # North waterfront, west to east
    (-122.4775, 37.8105),   # Fort Point
    (-122.4700, 37.8075),
    (-122.4650, 37.8060),   # Crissy Field
    (-122.4550, 37.8065),
    (-122.4430, 37.8072),
    (-122.4300, 37.8082),   # Fort Mason
    (-122.4235, 37.8078),   # Aquatic Park
    (-122.4180, 37.8090),
    (-122.4100, 37.8090),   # Pier 39
    (-122.4020, 37.8050),
    (-122.3960, 37.7990),
    (-122.3935, 37.7955),   # Ferry Building
    (-122.3880, 37.7860),
    (-122.3870, 37.7780),
    # inland / south filler and the outer (ocean) coast back north
    (-122.3800, 37.7000), (-122.5200, 37.7000),
    (-122.5115, 37.7600),   # Ocean Beach
    (-122.5060, 37.7850),   # Lands End
    (-122.4960, 37.7885),
    (-122.4830, 37.7930),   # Baker Beach
    (-122.4790, 37.8020),
]

EAST_BAY = [
    # Richmond down to Oakland, waterline west-facing
    (-122.3960, 37.9200),
    (-122.3900, 37.9120),   # Point Richmond
    (-122.3800, 37.9050),
    (-122.3700, 37.9010),
    (-122.3450, 37.8990),
    (-122.3300, 37.8985),
    (-122.3230, 37.8990),   # Point Isabel
    (-122.3195, 37.8900),
    (-122.3180, 37.8800),
    (-122.3175, 37.8710),   # Cesar Chavez Park / Berkeley Marina
    (-122.3155, 37.8655),   # marina entrance
    (-122.3140, 37.8600),
    (-122.3145, 37.8500),
    (-122.3130, 37.8410),   # Emeryville
    (-122.3160, 37.8330),
    (-122.3220, 37.8280),   # Bay Bridge toll plaza
    (-122.3250, 37.8150),
    (-122.3200, 37.8000),
    # inland filler
    (-122.2000, 37.7800), (-122.2000, 37.9200),
]

BROOKS_ISLAND = [
    (-122.3600, 37.8985), (-122.3500, 37.8990), (-122.3480, 37.8960),
    (-122.3580, 37.8955),
]

LAND_POLYGONS = {
    "marin": MARIN,
    "angel": ANGEL_ISLAND,
    "alcatraz": ALCATRAZ,
    "treasure": TREASURE_ISLAND,
    "yerba_buena": YERBA_BUENA,
    "san_francisco": SAN_FRANCISCO,
    "east_bay": EAST_BAY,
    "brooks": BROOKS_ISLAND,
}

# Effective height of each landmass in metres.  Used to scale how far downwind
# its wind shadow reaches -- Angel Island and the Marin headlands throw a real
# lee more than a mile long; flat Treasure Island barely throws one.
LAND_HEIGHT_M = {
    "marin": 280.0,          # Marin headlands / Wolfback Ridge
    "angel": 240.0,          # Mt Livermore
    "alcatraz": 40.0,
    "treasure": 8.0,
    "yerba_buena": 100.0,
    "san_francisco": 90.0,   # Presidio bluff / Fort Mason / city behind
    "east_bay": 15.0,        # low shoreline; hills are far inland
    "brooks": 50.0,
}

# --------------------------------------------------------------------------
# SI 10.0 OBSTRUCTIONS -- hard no-go.  Crossing one means retiring (SI 10.3).
#
# SI 10.1: "The lines between the closest associated shore to the following
# shall be considered obstructions."  Each entry is a polyline that starts at
# the charted object and runs to the shore.  If `shore` names a land polygon,
# the line is finished automatically at the nearest point of that polygon and
# pushed OBSTRUCTION_INLAND_M further inland so the barrier is sealed against
# the land mask at any grid resolution; if `shore` is None the polyline is
# used as written (it must already end inside land).
#
# Enforcement (geo.BayGrid + router.Router): every grid edge whose segment
# intersects one of these lines is removed from the graph, and every water
# cell within OBSTRUCTION_MARGIN_M of a line is made unnavigable.  That holds
# for all 16 move directions, so the barrier cannot be hopped diagonally.
# The same lines are exported to the browser app.
#
# Sources (all queried 11 Sep 2026):
#   ENC  = NOAA ENC Direct-to-GIS, harbour-scale cells US5OAKFG/FH/FF
#          (encdirect.noaa.gov), the digital form of charts 18649/18650/18653.
#   YRA  = YRA mark list, github.com/nkiesel/YRA-Marks (San_Francisco.gpx).
#   SSI  = StFYC / YRA sailing instructions text (Jessica Cup 2016 SI 7,
#          YRA Standing SI 10.1) -- descriptive only, no coordinates.
# --------------------------------------------------------------------------

OBSTRUCTION_MARGIN_M = 70.0     # safety strip either side of each line
OBSTRUCTION_INLAND_M = 200.0    # how far past the shoreline the line is sealed

OBSTRUCTIONS = [
    dict(
        key="berkeley_pier", si="10.1.a",
        name="Berkeley Pier ruins (daymark Fl R 4s 15ft 4M '2') to Berkeley shore",
        status="VERIFIED",
        source=("ENC US5OAKFH: light 37.847735N 122.360551W (COLOUR red, Fl 4s, "
                "4.6 m = 15 ft, 4 NM) at the outer end of the SLCONS 'Berkeley "
                "Pier' CONDTN=ruined polygons, which run to the standing fishing "
                "pier at 37.85952N 122.32702W and its root at 37.86287N "
                "122.31768W"),
        # daymark -> outer ruins -> standing pier -> root -> sealed inland
        line=[(-122.36055, 37.84774), (-122.33070, 37.85822),
              (-122.32702, 37.85952), (-122.31768, 37.86287),
              (-122.31300, 37.86450)],
        shore=None,
    ),
    dict(
        key="pt_blunt", si="10.1.b",
        name="Point Blunt buoy G '3' to Angel Island",
        status="VERIFIED",
        source=("ENC US5OAKFG lateral buoy 'San Francisco Bay North Channel "
                "Lighted Buoy 3', green pillar, Q G, at 37.850435N 122.417442W; "
                "YRA list agrees to 4 decimals"),
        line=[(-122.41744, 37.85043)],
        shore="angel",
    ),
    dict(
        key="alcatraz", si="10.1.c",
        name="Alcatraz Lighted Bell Buoy AZ (G/R/G, west end) to Alcatraz",
        status="VERIFIED",
        source=("ENC US5OAKFG lateral buoy 'Alcatraz Lighted Bell Buoy AZ', "
                "COLOUR 4,3,4 (green-red-green), at 37.827672N 122.428198W"),
        line=[(-122.42820, 37.82767)],
        shore="alcatraz",
    ),
    dict(
        key="h_beam", si="10.1.d",
        name="'H' beam piling off the old Water Quality Plant, W of StFYC, to shore",
        status="ESTIMATED (+-100 m) -- not a charted aid",
        source=("StFYC Jessica Cup SI 7: 'The H Beam piling (located "
                "approximately 200 yards west of the St. Francis YC)'.  ENC "
                "landmark 'stone tower' (StFYC) at 37.80756N 122.44385W; 200 yd "
                "west of the club's west end, just off East Beach, Crissy Field. "
                "Not in the ENC pile/obstruction/beacon layers.  Confirm by eye."),
        line=[(-122.4485, 37.8076)],
        shore="san_francisco",
    ),
    dict(
        key="anita_rock", si="10.1.e",
        name="Anita Rock, its light and the yellow 'AR' buoy, to Crissy Field",
        status="VERIFIED light; AR buoy position ESTIMATED from SI text",
        source=("ENC US5OAKFG special-purpose beacon 'Anita Rock Light' at "
                "37.808326N 122.453548W (Q W 20ft 4M, LNM 22/21), with awash "
                "rocks at 37.80821N 122.45405W / 37.80820N 122.45333W; YRA list "
                "'YRA-ANITA offset buoy' 37.80833N 122.45333W.  Jessica Cup SI: "
                "yellow buoy AR ~100 yd NW of the light -> (-122.4543, 37.8089)."),
        # AR buoy -> Anita Rock Light -> shore
        line=[(-122.45430, 37.80890), (-122.45355, 37.80833)],
        shore="san_francisco",
    ),
    dict(
        key="south_tower", si="10.1.f",
        name="Golden Gate Bridge South Tower to Fort Point (closes the inshore passage)",
        status="VERIFIED",
        source=("ENC US5OAKFF pylon area 'Golden Gate Bridge South Pier', "
                "polygon 37.81380-37.81422N 122.47737-122.47837W, centroid "
                "37.81401N 122.47787W; YRA list GGB-ST 37.813966N 122.478361W"),
        line=[(-122.47787, 37.81401)],
        shore="san_francisco",
    ),
]


# --------------------------------------------------------------------------
# SHIPPING LANES -- soft cost only (SI 9.1 / Inland Rule 9: do not impede).
#
# Regulated traffic lanes and precautionary areas from chart 18649/18650 as
# carried by the NOAA ENC approach cells US4CA2CJ / US4CA2BJ / US4CA11M
# (encdirect.noaa.gov, layer TSSLPT / PRCARE, queried 11 Sep 2026), hand-
# simplified to their corner points (the ENC rings hug Alcatraz and the
# waterfront vertex by vertex; that detail is below the grid resolution).
#
# These are NOT obstructions.  A sailboat may cross them; the router just
# charges `weight` x elapsed time while inside one, and flags the crossing so
# the tactician knows where a ship can force them off.  Weight 1.0 = flag
# only.  The Central Bay precautionary area covers essentially all of the
# water between Angel Island, Alcatraz and the Berkeley Circle buoy line, so
# it is flagged but not penalised -- there is nowhere else to sail.
# --------------------------------------------------------------------------

SHIPPING_ZONES = [
    dict(key="msc_wb", kind="lane", weight=1.10,
         name="Main Ship Channel westbound lane (33 CFR 167.405)",
         poly=[(-122.4979, 37.8175), (-122.5226, 37.8080), (-122.5294, 37.8052),
               (-122.5425, 37.8000), (-122.5705, 37.7888), (-122.5879, 37.7818),
               (-122.5864, 37.7790), (-122.5563, 37.7884), (-122.5295, 37.7974),
               (-122.5217, 37.8000), (-122.5171, 37.8015), (-122.4949, 37.8135)]),
    dict(key="msc_eb", kind="lane", weight=1.10,
         name="Main Ship Channel eastbound lane (33 CFR 167.405)",
         poly=[(-122.4949, 37.8135), (-122.5171, 37.8015), (-122.5217, 37.8000),
               (-122.5295, 37.7974), (-122.5563, 37.7884), (-122.5717, 37.7833),
               (-122.5701, 37.7802), (-122.5626, 37.7824), (-122.5295, 37.7920),
               (-122.5131, 37.7968), (-122.5069, 37.8000), (-122.4934, 37.8070),
               (-122.4909, 37.8083)]),
    dict(key="gg_pa", kind="precautionary", weight=1.05,
         name="Golden Gate precautionary area",
         poly=[(-122.4909, 37.8083), (-122.4979, 37.8175), (-122.4844, 37.8227),
               (-122.4793, 37.8246), (-122.4750, 37.8274), (-122.4698, 37.8313),
               (-122.4665, 37.8291), (-122.4629, 37.8265), (-122.4629, 37.8141),
               (-122.4782, 37.8142)]),
    dict(key="bay_wb", kind="lane", weight=1.10,
         name="SF Bay westbound traffic lane (north of Alcatraz)",
         poly=[(-122.4629, 37.8219), (-122.4628, 37.8261), (-122.4462, 37.8389),
               (-122.4392, 37.8398), (-122.4231, 37.8286), (-122.4256, 37.8284),
               (-122.4238, 37.8263)]),
    dict(key="bay_eb", kind="lane", weight=1.10,
         name="SF Bay eastbound traffic lane (south of Alcatraz)",
         poly=[(-122.4629, 37.8141), (-122.4629, 37.8213), (-122.4238, 37.8258),
               (-122.4204, 37.8257), (-122.4202, 37.8117), (-122.4397, 37.8137)]),
    dict(key="central_pa", kind="precautionary", weight=1.00,
         name="Central Bay precautionary area (flag only)",
         poly=[(-122.4202, 37.8120), (-122.4202, 37.8266), (-122.4231, 37.8286),
               (-122.4392, 37.8398), (-122.4406, 37.8462), (-122.4161, 37.8510),
               (-122.4161, 37.8642), (-122.3965, 37.8642), (-122.3965, 37.8228),
               (-122.3672, 37.8000), (-122.3594, 37.7839), (-122.3850, 37.7839),
               (-122.3950, 37.7950), (-122.4000, 37.8030), (-122.4100, 37.8080)]),
]

# Berkeley flats: shoal water inshore of roughly this line.  An Olson 25 draws
# ~4.5 ft, so this is a soft constraint -- treated as land only if
# `avoid_shoal=True`.
BERKELEY_SHOAL = [
    (-122.3175, 37.8800), (-122.3300, 37.8790), (-122.3350, 37.8700),
    (-122.3330, 37.8600), (-122.3200, 37.8520), (-122.3140, 37.8560),
    (-122.3175, 37.8700),
]


# --------------------------------------------------------------------------
# Point-in-polygon (vectorised ray casting) and mask construction
# --------------------------------------------------------------------------

def _points_in_poly(px, py, poly):
    """Vectorised even-odd point-in-polygon.  px, py are arrays in lon/lat."""
    poly = np.asarray(poly, dtype=float)
    xs, ys = poly[:, 0], poly[:, 1]
    inside = np.zeros(np.shape(px), dtype=bool)
    n = len(poly)
    j = n - 1
    for i in range(n):
        cond = ((ys[i] > py) != (ys[j] > py))
        # x coordinate of the edge at height py
        denom = np.where(ys[j] - ys[i] == 0, 1e-12, ys[j] - ys[i])
        xint = xs[i] + (py - ys[i]) * (xs[j] - xs[i]) / denom
        inside ^= cond & (px < xint)
        j = i
    return inside


# --------------------------------------------------------------------------
# Obstruction geometry helpers (exact, in projected metres -- not gridded)
# --------------------------------------------------------------------------

def _nearest_point_on_polygon(px, py, poly_xy):
    """Nearest point on the boundary of a closed ring to (px, py)."""
    P = np.asarray(poly_xy, float)
    A = P
    B = np.roll(P, -1, axis=0)
    V = B - A
    L2 = np.maximum((V ** 2).sum(axis=1), 1e-9)
    t = np.clip(((px - A[:, 0]) * V[:, 0] + (py - A[:, 1]) * V[:, 1]) / L2, 0, 1)
    Q = A + t[:, None] * V
    d2 = (Q[:, 0] - px) ** 2 + (Q[:, 1] - py) ** 2
    k = int(np.argmin(d2))
    return float(Q[k, 0]), float(Q[k, 1])


def obstruction_lines_xy(inland_m=OBSTRUCTION_INLAND_M):
    """
    The SI 10 obstruction lines as polylines in projected metres, each
    finished at the closest point ashore (nearest point of the named land
    polygon) and pushed `inland_m` further into the land so the barrier is
    sealed against the land mask.  Returns [(key, [(x, y), ...]), ...].
    """
    out = []
    for ob in OBSTRUCTIONS:
        pts = [tuple(map(float, ll_to_xy(lo, la))) for lo, la in ob["line"]]
        if ob["shore"] is not None:
            poly = [tuple(map(float, ll_to_xy(lo, la)))
                    for lo, la in LAND_POLYGONS[ob["shore"]]]
            qx, qy = _nearest_point_on_polygon(pts[-1][0], pts[-1][1], poly)
            dx, dy = qx - pts[-1][0], qy - pts[-1][1]
            L = max(np.hypot(dx, dy), 1e-9)
            pts.append((qx + dx / L * inland_m, qy + dy / L * inland_m))
        out.append((ob["key"], pts))
    return out


def obstruction_lines_ll(inland_m=OBSTRUCTION_INLAND_M):
    """Same as obstruction_lines_xy but in lon/lat (for export and plotting)."""
    out = []
    for key, pts in obstruction_lines_xy(inland_m):
        out.append((key, [tuple(map(float, xy_to_ll(x, y))) for x, y in pts]))
    return out


def _segments_intersect(ax, ay, bx, by, cx, cy, dx, dy):
    """
    Vectorised proper/improper segment intersection test: does segment AB
    intersect segment CD?  A, B may be arrays; C, D are scalars.  Touching
    counts as crossing (conservative -- the boat never gets the benefit of
    the doubt on an obstruction).
    """
    def orient(px, py, qx, qy, rx, ry):
        return (qx - px) * (ry - py) - (qy - py) * (rx - px)
    o1 = orient(ax, ay, bx, by, cx, cy)
    o2 = orient(ax, ay, bx, by, dx, dy)
    o3 = orient(cx, cy, dx, dy, ax, ay)
    o4 = orient(cx, cy, dx, dy, bx, by)
    return ((o1 * o2) <= 0) & ((o3 * o4) <= 0)


def _dist_to_segment(px, py, ax, ay, bx, by):
    vx, vy = bx - ax, by - ay
    L2 = max(vx * vx + vy * vy, 1e-9)
    t = np.clip(((px - ax) * vx + (py - ay) * vy) / L2, 0.0, 1.0)
    return np.hypot(px - (ax + t * vx), py - (ay + t * vy))


def track_crosses_obstruction(xy_points, lines_xy=None):
    """
    True if the polyline through xy_points (projected metres) crosses any
    obstruction line.  Used by the tests and by time_along_track.
    """
    lines_xy = lines_xy if lines_xy is not None else obstruction_lines_xy()
    P = np.asarray(xy_points, float)
    if len(P) < 2:
        return False
    ax, ay, bx, by = P[:-1, 0], P[:-1, 1], P[1:, 0], P[1:, 1]
    for _key, pts in lines_xy:
        for (cx, cy), (dx, dy) in zip(pts[:-1], pts[1:]):
            if np.any(_segments_intersect(ax, ay, bx, by, cx, cy, dx, dy)):
                return True
    return False


class BayGrid:
    """
    Regular grid in projected metres, with a water mask, the SI 10 obstruction
    barriers, and the shipping-lane zone index.

    Obstructions are enforced two ways, and both are needed:
      * `water` is False for any cell whose centre lies within
        `obstruction_margin_m` of an obstruction line, and
      * `blocked[r, c, j]` is True when the straight edge from cell (r, c) to
        its j-th neighbour (router.OFFSETS order) crosses an obstruction line.
    The first alone leaks -- a 250 m grid with 2:1 knight moves can hop a
    70 m strip -- and the second alone would let a route graze the line.
    """

    def __init__(self, res_m=250.0, avoid_shoal=False, shore_buffer_m=60.0,
                 obstructions=True, obstruction_margin_m=OBSTRUCTION_MARGIN_M,
                 include_pier=None):
        # `include_pier` is accepted for backward compatibility; the pier is
        # now SI obstruction 10.1.a and is controlled by `obstructions`.
        if include_pier is not None:
            obstructions = bool(include_pier)
        self.res = float(res_m)
        d = DOMAIN
        x0, y0 = ll_to_xy(d["lon_min"], d["lat_min"])
        x1, y1 = ll_to_xy(d["lon_max"], d["lat_max"])
        self.nx = int(np.ceil((x1 - x0) / self.res))
        self.ny = int(np.ceil((y1 - y0) / self.res))
        self.x = x0 + (np.arange(self.nx) + 0.5) * self.res
        self.y = y0 + (np.arange(self.ny) + 0.5) * self.res
        self.X, self.Y = np.meshgrid(self.x, self.y)          # (ny, nx)
        self.LON, self.LAT = xy_to_ll(self.X, self.Y)

        land = np.zeros((self.ny, self.nx), dtype=bool)
        height = np.zeros((self.ny, self.nx), dtype=float)
        for name, poly in LAND_POLYGONS.items():
            m = _points_in_poly(self.LON, self.LAT, poly)
            land |= m
            height = np.maximum(height, m * LAND_HEIGHT_M.get(name, 50.0))
        if avoid_shoal:
            land |= _points_in_poly(self.LON, self.LAT, BERKELEY_SHOAL)

        self.land = land
        self.height = height
        self.water = ~land

        # Distance (in metres) from each water cell to the nearest land cell.
        # Used for shoreline current relief and shoreline wind shadowing.
        from scipy import ndimage
        self.dist_to_land = ndimage.distance_transform_edt(
            self.water, sampling=self.res)

        # Shave a thin strip off the waterline so routes do not graze rocks.
        # NOTE: this narrows `water` (navigability) only.  `land` stays the
        # polygon land, because the wind lee ray-march stops at the first land
        # cell and reads its height -- folding a zero-height buffer strip into
        # `land` would stop every ray one cell early and erase the shadows.
        if shore_buffer_m > 0:
            self.water = self.water & (self.dist_to_land >= shore_buffer_m)

        # -- SI 10 obstructions ------------------------------------------
        from .router import OFFSETS as _OFFSETS
        self.obstruction_lines = obstruction_lines_xy() if obstructions else []
        self.obstruction_margin = float(obstruction_margin_m)
        self.obstructed = np.zeros((self.ny, self.nx), dtype=bool)
        self.blocked = np.zeros((self.ny, self.nx, len(_OFFSETS)), dtype=bool)
        for _key, pts in self.obstruction_lines:
            for (ax, ay), (bx, by) in zip(pts[:-1], pts[1:]):
                self.obstructed |= (_dist_to_segment(self.X, self.Y, ax, ay, bx, by)
                                    < self.obstruction_margin)
                for j, (dc, dr) in enumerate(_OFFSETS):
                    ex = self.X + dc * self.res
                    ey = self.Y + dr * self.res
                    self.blocked[:, :, j] |= _segments_intersect(
                        self.X, self.Y, ex, ey, ax, ay, bx, by)
        self.water = self.water & ~self.obstructed

        # -- shipping zones (soft) ----------------------------------------
        self.zone = np.full((self.ny, self.nx), -1, dtype=np.int16)
        self.zone_weight = np.ones(len(SHIPPING_ZONES) + 1)
        for k, z in enumerate(SHIPPING_ZONES):
            m = _points_in_poly(self.LON, self.LAT, z["poly"]) & (self.zone < 0)
            self.zone[m] = k
            self.zone_weight[k] = float(z["weight"])
        self.zone_weight[-1] = 1.0                 # index -1 -> no zone

    # -- index helpers ----------------------------------------------------

    def nearest_water_cell(self, lon, lat):
        """(row, col) of the nearest water cell to a lon/lat."""
        x, y = ll_to_xy(lon, lat)
        d2 = (self.X - x) ** 2 + (self.Y - y) ** 2
        d2 = np.where(self.water, d2, np.inf)
        idx = int(np.argmin(d2))
        return np.unravel_index(idx, d2.shape)

    def cell_center_ll(self, rc):
        r, c = rc
        return float(self.LON[r, c]), float(self.LAT[r, c])

    def cell_xy(self, rc):
        r, c = rc
        return float(self.X[r, c]), float(self.Y[r, c])

    @property
    def extent_ll(self):
        return (float(self.LON.min()), float(self.LON.max()),
                float(self.LAT.min()), float(self.LAT.max()))


def water_distance_field(grid, lon, lat):
    """
    Geodesic distance *through water* from a lon/lat seed to every water cell,
    in metres.  Used so that current/wind interpolation does not leak across
    headlands (e.g. Sausalito should not borrow from the Cityfront).
    """
    import heapq
    ny, nx = grid.ny, grid.nx
    D = np.full((ny, nx), np.inf)
    r0, c0 = grid.nearest_water_cell(lon, lat)
    D[r0, c0] = 0.0
    heap = [(0.0, r0 * nx + c0)]
    nbrs = [(-1, 0, 1.0), (1, 0, 1.0), (0, -1, 1.0), (0, 1, 1.0),
            (-1, -1, 1.4142), (-1, 1, 1.4142), (1, -1, 1.4142), (1, 1, 1.4142)]
    water = grid.water
    res = grid.res
    while heap:
        d, k = heapq.heappop(heap)
        r, c = divmod(k, nx)
        if d > D[r, c]:
            continue
        for dr, dc, w in nbrs:
            rr, cc = r + dr, c + dc
            if 0 <= rr < ny and 0 <= cc < nx and water[rr, cc]:
                nd = d + w * res
                if nd < D[rr, cc]:
                    D[rr, cc] = nd
                    heapq.heappush(heap, (nd, rr * nx + cc))
    return D

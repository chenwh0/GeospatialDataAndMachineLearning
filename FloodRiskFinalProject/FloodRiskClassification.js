 // 1. SETUP
var countries = ee.FeatureCollection("FAO/GAUL/2015/level0");
var taiwan = countries.filter(ee.Filter.eq("ADM0_NAME", "Taiwan"));
var REGION = taiwan.geometry();

var startDate = "2025-05-01";
var endDate   = "2025-12-01";

// Common classification grid (m)
var GRID_SCALE = 150;         
var GRID_CRS   = "EPSG:4326";
var GRID_PROJ  = ee.Projection(GRID_CRS).atScale(GRID_SCALE);

// 2. HELPER FUNCTIONS
// Continuous to mean on grid
function toGridMean(img) {
  return img
    .setDefaultProjection(GRID_PROJ)
    .reduceResolution({
      reducer: ee.Reducer.mean(),
      maxPixels: 1024
    })
    .reproject(GRID_PROJ);
}

// Categorical to mode on grid
function toGridMode(img) {
  return img
    .setDefaultProjection(GRID_PROJ)
    .reduceResolution({
      reducer: ee.Reducer.mode(),
      maxPixels: 1024
    })
    .reproject(GRID_PROJ);
}

// Extensive → sum on grid
function toGridSum(img) {
  return img
    .setDefaultProjection(GRID_PROJ)
    .reduceResolution({
      reducer: ee.Reducer.sum(),
      maxPixels: 1024
    })
    .reproject(GRID_PROJ);
}

// Min–max normalization (0–1) at grid scale
function normalize(img) {
  var band = img.bandNames().get(0);
  var stats = img.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: REGION,
    scale: GRID_SCALE,
    maxPixels: 1e8,
    bestEffort: true
  });
  var min = ee.Number(stats.get(ee.String(band).cat("_min")));
  var max = ee.Number(stats.get(ee.String(band).cat("_max")));
  return img.subtract(min).divide(max.subtract(min)).clamp(0, 1);
}

// PREPARE FEATURES INPUT DATA

// HAZARD INDICES
// Precipitation (CHIRPS, accumulated)
var precipitation = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
  .filterDate(startDate, endDate)
  .filterBounds(REGION)
  .sum();
precipitation = toGridMean(precipitation).rename("precipitation");

// elevation & slope
var elevation = toGridMean(ee.Image("USGS/SRTMGL1_003")).rename("elevation");
var slope = ee.Terrain.slope(elevation).rename("slope");
// TWI proxy
var slopeRad = slope.multiply(Math.PI).divide(180);
var upslope = elevation.reduceNeighborhood({
  reducer: ee.Reducer.sum(),
  kernel: ee.Kernel.square(3)
});
var twi = upslope
  .add(1)
  .divide(slopeRad.tan().add(0.001))
  .log()
  .rename("twi");
twi = toGridMean(twi);

// Runoff (FLDAS Storm runoff 95th percentile)
var storm_runoff_p95 = ee.ImageCollection("NASA/FLDAS/NOAH01/C/GL/M/V001")
  .filterDate(startDate, endDate)
  .select("Qs_tavg")
  .reduce(ee.Reducer.percentile([95]))
  .rename("storm_runoff_p95");
storm_runoff_p95 = toGridMean(storm_runoff_p95);


// EXPOSURE INDICES
// Land cover (ESA WorldCover)
var landcover = toGridMode(
  ee.Image("ESA/WorldCover/v200/2021").select("Map")
).rename("landcover");

var urban    = landcover.eq(50).rename("urban");
var cropland = landcover.eq(40).rename("cropland");

// Population (WorldPop)
var population = ee.ImageCollection("WorldPop/GP/100m/pop")
  .filterBounds(REGION)
  .filter(ee.Filter.eq("year", 2020))
  .select("population")
  .mosaic();
population = toGridSum(population).rename("population");


// VULNERABILITY INDICES
// SMAP wetness
var smap = ee.ImageCollection("NASA/SMAP/SPL4SMGP/008")
  .filterBounds(REGION)
  .first();
var rootzone = toGridMean(smap.select("sm_rootzone_wetness")).rename("rootzone");
var profile = toGridMean(smap.select("sm_profile_wetness")).rename("profile");

// 6. NORMALIZATION (0–1)
var precipitation_n = normalize(precipitation);
var elevation_n = normalize(ee.Image(1).subtract(elevation));
var slope_n = normalize(ee.Image(1).subtract(slope));
var twi_n = normalize(twi);
var storm_runoff_n = normalize(storm_runoff_p95);
var population_n = normalize(population);
var root_n = normalize(rootzone);
var profile_n = normalize(profile);

// 7. FEATURE STACK
var features = ee.Image.cat([
  slope_n.rename("slope"),
  elevation_n.rename("elevation"),
  twi_n.rename("twi"),
  precipitation_n.rename("precipitation"),
  storm_runoff_n.rename("runoff"),
  population_n.rename("population"),
  urban.multiply(1).rename("urban"),
  cropland.multiply(0.6).rename("cropland"),
  root_n.rename("rootzone"),
  profile_n.rename("profile")
]);

// 8. FLOOD RISK INDEX
var flood_risk =
    slope_n.multiply(0.13)
  .add(elevation_n.multiply(0.13))
  .add(twi_n.multiply(0.15))
  .add(precipitation_n.multiply(0.20))
  .add(storm_runoff_n.multiply(0.15))
  .add(population_n.multiply(0.10))
  .add(urban.multiply(0.04))
  .add(cropland.multiply(0.02))
  .add(root_n.multiply(0.04))
  .add(profile_n.multiply(0.04))
  .rename("FloodRiskIndex");


// 9. CLASSIFICATION (PERCENTILES)
var pct = flood_risk.reduceRegion({
  reducer: ee.Reducer.percentile([20, 40, 60, 80]),
  geometry: REGION,
  scale: GRID_SCALE,
  maxPixels: 1e8,
  bestEffort: true
});

var p20 = ee.Number(pct.get("FloodRiskIndex_p20"));
var p40 = ee.Number(pct.get("FloodRiskIndex_p40"));
var p60 = ee.Number(pct.get("FloodRiskIndex_p60"));
var p80 = ee.Number(pct.get("FloodRiskIndex_p80"));

var predicted_risk = ee.Image(1)
  .where(flood_risk.gte(p20), 2)
  .where(flood_risk.gte(p40), 3)
  .where(flood_risk.gte(p60), 4)
  .where(flood_risk.gte(p80), 5)
  .rename("FloodRiskClass");

// 10. PREPARE TRUE DATA
// Get 100-year flood reocurrence - JRC, 90m resolution (true data)
var true_dataset = ee.ImageCollection("JRC/CEMS_GLOFAS/FloodHazard/v2_1");
var true_risk_RP100_depth = true_dataset.select("RP100_depth")
  .filterBounds(REGION)
  .mosaic(); // JRC only has data for March 16, 2024
true_risk_RP100_depth = toGridMean(true_risk_RP100_depth);

var true_pct = true_risk_RP100_depth.reduceRegion({
  reducer: ee.Reducer.percentile([20, 40, 60, 80]),
  geometry: REGION,
  scale: GRID_SCALE,
  maxPixels: 1e8,
  bestEffort: true
});

var tp20 = ee.Number(true_pct.get("RP100_depth_p20"));
var tp40 = ee.Number(true_pct.get("RP100_depth_p40"));
var tp60 = ee.Number(true_pct.get("RP100_depth_p60"));
var tp80 = ee.Number(true_pct.get("RP100_depth_p80"));

var true_risk = ee.Image(1)
  .where(true_risk_RP100_depth.gte(tp20), 2)
  .where(true_risk_RP100_depth.gte(tp40), 3)
  .where(true_risk_RP100_depth.gte(tp60), 4)
  .where(true_risk_RP100_depth.gte(tp80), 5)
  .rename("TrueFloodRiskClass");
  
  
// COMPARE USING CONFUSION MATRIX
// Stack predicted vs true
var stack = predicted_risk.addBands(true_risk);

// Generate confusion matrix
var sample = stack.sample({
  region: REGION,
  scale: GRID_SCALE,
  numPixels: 1e6
});

var confusion_matrix = sample.errorMatrix("TrueFloodRiskClass", "FloodRiskClass");

print('Confusion Matrix:', confusion_matrix);
print('Overall Accuracy:', confusion_matrix.accuracy());
print('Kappa:', confusion_matrix.kappa());
print('Producers Accuracy:', confusion_matrix.producersAccuracy());
print('Users Accuracy:', confusion_matrix.consumersAccuracy());

// VISUALIZATION 
// Color palette & mapping
var colors = [
  "#f7fbff", // 1 Very Low
  "#c6dbef", // 2 Low
  "#6baed6", // 3 Moderate
  "#2171b5", // 4 High
  "#08306b"  // 5 Very High
];

var risk_map_colors = {
  min: 1,
  max: 5,
  palette: colors
};

// Map Legend
var legend = ui.Panel({style: { position: "bottom-right", padding: "8px" }});
legend.add(ui.Label({value: "Flood Risk Classes", style: { fontWeight: "bold" }}));

var names = [
  "1 Very Low",
  "2 Low",
  "3 Moderate",
  "4 High",
  "5 Very High"
];

for (var i = 0; i < names.length; i++) {
  var row = ui.Panel({layout: ui.Panel.Layout.flow("horizontal")});
  var colorBox = ui.Label("", {backgroundColor: colors[i], padding: "8px", margin: "0 6px 0 0"});
  var label = ui.Label(names[i]);
  row.add(colorBox);
  row.add(label);
  legend.add(row);
}

// Bar chart: Area by Class
var areaImage = ee.Image.pixelArea()
  .divide(1e6)          // m² → km²
  .addBands(predicted_risk.clip(REGION)); // clip to Taiwan

var areaChart = ui.Chart.image.byClass({
  image: areaImage,
  classBand: "FloodRiskClass",
  region: REGION,       // restrict chart to Taiwan
  scale: GRID_SCALE,
  reducer: ee.Reducer.sum()
})
.setChartType("ColumnChart")
.setOptions({
  title: "Flood Risk Area by Class (km²)",
  hAxis: { title: "Flood Risk Class" },
  vAxis: { title: "Area (km²)" },
  legend: { position: "none" },
  colors: colors
});

print(areaChart);

// Map
var white_background = ee.Image.constant(1).visualize({palette: ["white"]});
Map.addLayer(white_background, {}, "White background");
Map.addLayer(REGION, {color: "black"}, "study region", false);
Map.centerObject(REGION, 8);
Map.add(legend);
Map.addLayer(predicted_risk.clip(REGION), risk_map_colors, "Flood Risk Classes");

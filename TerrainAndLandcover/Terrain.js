var srtm = ee.Image("USGS/SRTMGL1_003");  // 30 m SRTM
var region = ee.Geometry.Rectangle([-120.65,35.65, -118.014,37.48]); // Sierra Nevada and Central Valley transect

// Get srtm of region
var srtm_region = srtm.clip(region);


// Task 6.4 - Flood risk mapping
// Get smoothed srtm of region
var region_smoothed = srtm_region.focal_mean({radius: 1, units: "pixels"})
               .focal_min({radius: 1, units: "pixels"});
// Get slope
var slopeDeg = ee.Terrain.slope(region_smoothed);
var slopeRad = slopeDeg.multiply(Math.PI).divide(180);

// 3) Proxy flow accumulation (demo-scale)
// Note: Earth Engine doesn't natively expose full D8/D-infinity at global scale. This proxy uses a coarse neighborhood sum as an upslope "area" indicator.
// For production hydrology, export the DEM and compute flow accumulation in SAGA/QGIS/Arc.
var kernel = ee.Kernel.square({radius: 3, units: "pixels"});  // ~210 m radius at 30 m
var upslopeProxy = region_smoothed.reduceNeighborhood({reducer: ee.Reducer.sum(), kernel: kernel});

// 4) Compute TWI (stabilized)
// TWI = ln( a / tan(beta) ), add small epsilon to avoid division by zero
var eps = 0.001;
var twi = upslopeProxy.add(1).divide(slopeRad.tan().add(eps)).log().rename("TWI");



// Flood-risk conditions
var low_slope_threshold = 5; // Lower slope (degrees) increases flood risk
var low_elevation_threshold = 40; // Lower elevation (meters) increases river-flood risk
var high_TWI_threshold = 12; // Higher TWI can increase flood risk

var low_slope_mask = slopeDeg.lt(low_slope_threshold).rename("lowSlope");
var low_elevation_mask = srtm_region.lt(low_elevation_threshold).rename("lowElevation");
var high_TWI_mask = twi.gt(high_TWI_threshold).rename("highTWI");

var slope_risk = low_slope_mask.multiply(1);
var elevation_risk = low_elevation_mask.multiply(1);
var twi_risk = high_TWI_mask.multiply(1);

var total_risk = slope_risk.add(elevation_risk).add(twi_risk).rename("RiskClass"); // 0-3 range. 0 = no risk criterias met & 3 = all risk criterias met



// Create a histogram chart for total_risk
var counts = total_risk.reduceRegion({
  reducer: ee.Reducer.frequencyHistogram(), // counts pixels per value
  geometry: region,
  scale: 30,                   // DEM resolution
  maxPixels: 1e8                // increase limit
});
print('Pixel counts', counts);

// Create a histogram chart
var chart = ui.Chart.image.histogram({
  image: total_risk,
  region: region,  // your study area geometry
  scale: 30,       // match your DEM resolution
  maxPixels: 1e8,
  maxBuckets: 4    // risk levels 0,1,2,3
})
.setOptions({
  title: 'Total Flood Risk Distribution',
  hAxis: {title: 'Risk Level (0=Low, 3=High)'},
  vAxis: {title: 'Number of Pixels'},
  colors: ['#ff0000'],
  legend: {position: 'none'}
});
print(chart);
/*

// Visualize
var risk_visualize_params = {
  min: 0, 
  max: 3, 
  palette: [
    "white", // 0, = no/low flood risk
    "caf0f8", 
    "0077b6", 
    "03045e" // 3, darkblue = high flood risk
    ]
};

Map.addLayer(region, {}, "study region", false);
Map.centerObject(region, 8);
Map.addLayer(twi.updateMask(twi), {min: 4, max: 12, palette: ['red','yellow','lightgreen','darkgreen','blue']}, 'TWI (proxy)');
Map.addLayer(total_risk, risk_visualize_params, "Flood risk");

/*
// Task 6.3 - Integrate land cover class with slope
// Get slope
var terrain = ee.Algorithms.Terrain(srtm_region);
var slope = terrain.select("slope"); // rate of elevation change (steepness of an incline)
// Get cropland classes
var nlcd2021 = ee.Image("USGS/NLCD_RELEASES/2021_REL/NLCD/2021").clip(region);
var landcover = nlcd2021.select("landcover");
var nlcdClasses = [11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95];

// Group by landcover class
var reducer = ee.Reducer.sum()
  .combine(ee.Reducer.count(), "", true)   // sharedInputs = true
  .group({
    groupField: 1,    // landcover band index
    groupName: "class"
  });
var stats = slope.addBands(landcover).reduceRegion({
  reducer: reducer,
  geometry: region,
  scale: 30,
  maxPixels: 1e13
});
var groups = ee.List(stats.get('groups')); 


// Compute each group's mean = sum / count
var table = ee.FeatureCollection(
  groups.map(function(item) {
    var d = ee.Dictionary(item);
    var sum = ee.Number(d.get('sum'));
    var count = ee.Number(d.get('count'));
    var mean = ee.Algorithms.If(count.gt(0), sum.divide(count), null);
    return ee.Feature(null, {
      NLCD_Class: d.get('class'),
      Slope_Mean: mean
    });
  })
);

// Remove any classes with zero pixels
var table_clean = table.filter(ee.Filter.notNull(['Slope_Mean']));

// Isolate slope of only cropland
var isCropland = landcover.eq(82);
var slope_only_cropland = slope.updateMask(isCropland);

// Visualize
print("Mean slope per NLCD class", table_clean);
var chart = ui.Chart.feature.byFeature({
  features: table_clean.sort('NLCD_Class'),
  xProperty: 'NLCD_Class',
  yProperties: ['Slope_Mean']
})
.setChartType('ColumnChart')
.setOptions({
  title: 'Mean Slope by NLCD Class',
  hAxis: {title: 'NLCD Class'},
  vAxis: {title: 'Slope (degrees)'},
  legend: {position: 'none'}
});
print(chart);

var nlcdPalette = [
  '476ba1', // 11 Open Water
  'd1def8', // 12 Perennial Ice/Snow  
  'dec5c5', // 21 Developed, Open Space
  'd99282', // 22 Developed, Low Intensity
  'eb0000', // 23 Developed, Medium Intensity
  'ab0000', // 24 Developed, High Intensity
  'b3ac9f', // 31 Barren
  '68ab5f', // 41 Deciduous Forest
  '1c5f2c', // 42 Evergreen Forest
  'b5c58f', // 43 Mixed Forest
  'af963c', // 52 Shrub/Scrub
  'ccb879', // 71 Grassland/Herbaceous
  'dfdfc2', // 81 Pasture/Hay
  'd1d182', // 82 Cultivated Crops
  'a3cc51', // 90 Woody Wetlands
  '82b3d1'  // 95 Emergent Herbaceous Wetlands
];

Map.addLayer(region, {}, "study region", false);
Map.centerObject(region, 8);
Map.addLayer(landcover, {min: 0, max: 95, palette: nlcdPalette}, 'NLCD 2021');
Map.addLayer(isCropland.selfMask(), {palette: ['#ffd700']}, 'Cropland (NLCD 82)', false);

// Visualize slope of only cropland
var slopeVisualizeParams = {
  min: 0,            
  max: 45,
  palette: ["green","yellow","red"] 
};

Map.addLayer(region);
Map.centerObject(region, 8);
Map.addLayer(slope_only_cropland, slopeVisualizeParams, "Slope of Cropland");
*/

/*
// Task 6.2
// Calculate variables based on terrain data - slope, aspect, hillshade
var terrain = ee.Algorithms.Terrain(srtm_region);
var slope = terrain.select("slope"); // rate of elevation change (steepness of an incline)
var aspect = terrain.select("aspect"); // Compass direction of max slope (direction water will flow)
var hillshade = ee.Terrain.hillshade(srtm_region, 315, 45); // illumination of the terrain by a hypothetical sun. Good for black-white visualization

function getStats(image, name) {
  var stats = image.reduceRegion({
    reducer: ee.Reducer.mean()
      .combine(ee.Reducer.minMax(), "", true)
      .combine(ee.Reducer.stdDev(), "", true),
    geometry: region,
    scale: 30,
    maxPixels: 1e13
  });
  
  print(name, "stats:", stats);
}

getStats(slope, "slope");
getStats(aspect, "aspect");
getStats(hillshade, "hillshade");

// Visualize
var slopeVisualizeParams = {
  min: 0,            
  max: 45,
  palette: ["green","yellow","red"] 
};

var aspectVisualizeParams = {
  min: 0,            
  max: 360,
  palette: ["blue","green","yellow", "orange", "red", "purple"] // blue = north, green = east, red = south, purple = west
};

var hillshadeVisualizeParams = {min: 100, max: 255};

Map.addLayer(region, {}, "study region", false);
Map.centerObject(region, 8);
Map.addLayer(slope, slopeVisualizeParams, "Slope (degrees)");
Map.addLayer(aspect, aspectVisualizeParams, "Aspect (degrees)");
Map.addLayer(hillshade, hillshadeVisualizeParams, "Hillshade");
*/

/*
// Task 6.1
// Fetch SRTM
var min_max_stats = srtm_region.reduceRegion({
  reducer: ee.Reducer.minMax(),
  geometry: region,
  scale: 30,
  maxPixels: 1e13
});

print("Min & max elevation stats for Sierra Nevada and Central Valley region (m):", min_max_stats);

// Visualize
var elevationVisualizeParams = {
  min: 0,            // use range given by min_max_stats
  max: 4400,
  palette: ["2166ac","5aae61","a6d96a","ffffbf","fdae61","f46d43","d73027"] // blue, green, light green, yellow, orange, red-orange, red
};

Map.addLayer(region, {}, "study region", false);
Map.centerObject(region, null);
Map.addLayer(srtm_region, elevationVisualizeParams, "SRTM study region");

// Export geoTIFF to Google Drive
Export.image.toDrive({
  image: srtm_region,
  description: "SRTM_30m_CentralValley",
  fileNamePrefix: "SRTM_30m_CentralValley",
  region: region,
  scale: 30,
  crs: "EPSG:4326",
  fileFormat: "GeoTIFF",
  maxPixels: 1e13,
  formatOptions: {cloudOptimized: true}
});
*/
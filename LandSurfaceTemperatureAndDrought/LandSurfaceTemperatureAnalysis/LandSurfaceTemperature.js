/*
var coordinates = ee.Geometry.Point([-118.243683, 34.052235]); // Los Angeles, California
var states = ee.FeatureCollection("TIGER/2018/States");
var california = states.filter(ee.Filter.eq("NAME", "California"));

// Get MODIS Terra LST 8-day composites for 2023
var modis = ee.ImageCollection("MODIS/061/MOD11A2")
  .filterDate("2023-01-01", "2024-01-01")
  .select(["LST_Day_1km", "LST_Night_1km"]);

// Convert to Celsius & apply quality control
var convertToCelsius = function(img) {
  var lstDay = img.select("LST_Day_1km").multiply(0.02).subtract(273.15);
  var lstNight = img.select("LST_Night_1km").multiply(0.02).subtract(273.15);
  
  return img.addBands(lstDay.rename("LST_Day_C"))
    .addBands(lstNight.rename("LST_Night_C"))
    .copyProperties(img, ["system:time_start"]);
};

var modisCelsius = modis.map(convertToCelsius)

var lstTimeSeries = modisCelsius.map(function(img) {
  var avgDay = img.select("LST_Day_C").reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: coordinates,
    scale: 1000,
    maxPixels: 1e9
  }).get("LST_Day_C");
  var avgNight = img.select("LST_Night_C").reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: coordinates,
    scale: 1000,
    maxPixels: 1e9
  }).get("LST_Night_C");
  
  return ee.Feature(null, {
    date: ee.Date(img.get("system:time_start")).format("YYYY-MM-dd"),
    day_temp: avgDay,
    night_temp: avgNight,
    date_milliseconds: img.get("system:time_start")
  });
});

// Create feature collections for day & night
var dayCollection = lstTimeSeries.map(function(feat) {
  return feat.set("temperature", feat.get("day_temp"))
             .set("series", "Day LST");
});

var nightCollection = lstTimeSeries.map(function(feat) {
  return feat.set("temperature", feat.get("night_temp"))
             .set("series", "Night LST");
});

var combinedCollection = dayCollection.merge(nightCollection);

// Create dual-series time series chart
var chart = ui.Chart.feature.groups(
  combinedCollection, "date", "temperature", "series")
  .setChartType("LineChart")
  .setOptions({
    title: "MODIS LST Seasonal Cycle - Los Angeles, California - 2023",
    hAxis: {
      title: "Date",
      titleTextStyle: {italic: false, bold: true}
    },
    vAxis: {
      title: "Land Surface Temperature (°C)",
      titleTextStyle: {italic: false, bold: true}
    },
    lineWidth: 2,
    series: {
      0: {color: "#d62728"}, // Day LST (red)
      1: {color: "#1f77b4"}  // Night LST (blue)
    },
    legend: {position: "top", alignment: "center"}
});

print(chart);

// Seasonal statistics
var seasonal = lstTimeSeries.aggregate_array("day_temp");
print("2023 Day LST seasonal statistics", seasonal.reduce(ee.Reducer.minMax()));
var annualAvgDay = modisCelsius.select("LST_Day_C").mean();
var annualAvgNight = modisCelsius.select("LST_Night_C").mean();
var visualizationVars = {
  min: 0, max: 35,
  palette: ["#000080", "#0000FF", "#0080FF", "#00FFFF", "#80FF80", 
            "#FFFF00", "#FF8000", "#FF0000", "#800000"]
};

Map.centerObject(california, 7);
Map.addLayer(annualAvgDay.clip(california), visualizationVars, "MODIS Annual Mean Day LST 2023 (°C)");
Map.addLayer(annualAvgNight.clip(california), visualizationVars, "MODIS Annual Mean Night LST 2023 (°C)");
Map.addLayer(coordinates, {color: "black", pointRadius: 5}, "Los Angeles, California");
*/

// ========================================================================
// URBAN HEAT ISLAND ANALYSIS - LANDSAT & SENTINEL-2
// ========================================================================

// Define Los Angeles metropolitan area

var losAngeles = ee.Geometry.Rectangle([-118.7, 33.7, -117.6, 34.5]);

// ========================================================================
// UTILITY FUNCTIONS
// ========================================================================

// Function to calculate Landsat LST
var calculateLandsatLST = function(image) {
  // Apply Collection 2 scaling to optical bands
  var optical = image.select(["SR_B2","SR_B3","SR_B4","SR_B5","SR_B6"])
                     .multiply(0.0000275).add(-0.2);

  // Scale thermal band to Kelvin
  var thermal = image.select("ST_B10").multiply(0.00341802).add(149.0);

  // Calculate NDVI for vegetation analysis
  var ndvi = optical.normalizedDifference(["SR_B5", "SR_B4"]).rename("NDVI");

  // Calculate proportion of vegetation for emissivity
  var pv = ndvi.subtract(0.2).divide(0.6).pow(2).clamp(0, 1);

  // Surface emissivity calculation
  var emissivity = pv.multiply(0.004).add(0.986).rename("Emissivity");

  // Single-channel LST algorithm
  var lambda = 11.5e-6; // Band 10 wavelength (meters)
  var rho = 1.438e-2;   // h*c/σ constant

  var lst = thermal.expression(
    "TB / (1 + (lambda * TB / rho) * log(eps))", {
      "TB": thermal,
      "lambda": lambda,
      "rho": rho,
      "eps": emissivity
    }
  ).subtract(273.15).rename("LST"); // Convert to Celsius

  return image.addBands([lst, ndvi, emissivity])
              .copyProperties(image, ["system:time_start", "CLOUD_COVER"]);
};

// cloud masking for Landsat
var maskLandsatClouds = function(image) {
  var qa = image.select("QA_PIXEL");
  var cloud = qa.bitwiseAnd(1 << 3).neq(0);       // Cloud
  var cloudShadow = qa.bitwiseAnd(1 << 4).neq(0); // Cloud shadow
  var snow = qa.bitwiseAnd(1 << 5).neq(0);        // Snow

  var mask = cloud.or(cloudShadow).or(snow).not();
  return image.updateMask(mask);
};

// Sentinel-2 preprocessing
var processSentinel2 = function(image) {
  // Scale Sentinel-2 reflectance values
  var optical = image.select(["B2","B3","B4","B8","B11"])
                     .multiply(0.0001);

  // Calculate NDVI
  var ndvi = optical.normalizedDifference(["B8", "B4"]).rename("NDVI");

  // LST estimation using NDVI and SWIR relationship
  // Empirical model: LST ≈ base_temp + SWIR_effect - NDVI_cooling
  var swir1 = optical.select("B11");
  var lstEstimate = swir1.multiply(60)      // SWIR warming effect
                         .subtract(ndvi.multiply(12))  // NDVI cooling effect
                         .add(30)          // Base temperature
                         .rename("Estimated_LST");

  // cloud masking
  var cloudProb = image.select("MSK_CLDPRB");
  var cloudMask = cloudProb.lt(30); // 30% cloud probability threshold

  return image.addBands([ndvi, lstEstimate])
              .updateMask(cloudMask)
              .copyProperties(image, ["system:time_start"]);
};

// UHI calculation using temperature percentiles
var calculateUHI = function(lstImage, geometry, bandName) {
  // Calculate basic statistics
  var stats = lstImage.select(bandName).reduceRegion({
    reducer: ee.Reducer.percentile([10, 90]).combine({
      reducer2: ee.Reducer.mean(),
      sharedInputs: true
    }),
    geometry: geometry,
    scale: 100,
    maxPixels: 1e9,
    bestEffort: true
  });

  var meanTemp = ee.Number(stats.get(bandName + "_mean"));
  var hotThreshold = ee.Number(stats.get(bandName + "_p90"));
  var coolThreshold = ee.Number(stats.get(bandName + "_p10"));

  // Calculate UHI intensity (hot areas - cool areas)
  var hotAreas = lstImage.select(bandName).updateMask(
    lstImage.select(bandName).gt(hotThreshold));
  var coolAreas = lstImage.select(bandName).updateMask(
    lstImage.select(bandName).lt(coolThreshold));

  var hotMean = hotAreas.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geometry,
    scale: 100,
    maxPixels: 1e9
  });

  var coolMean = coolAreas.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geometry,
    scale: 100,
    maxPixels: 1e9
  });

  var uhiIntensity = ee.Number(hotMean.get(bandName))
                      .subtract(ee.Number(coolMean.get(bandName)));

  // Create UHI map (temperature relative to mean)
  var uhiMap = lstImage.select(bandName).subtract(meanTemp)
                      .rename(bandName + "_UHI");

  return {
    uhiMap: uhiMap,
    uhiIntensity: uhiIntensity,
    meanTemp: meanTemp,
    hotTemp: hotMean.get(bandName),
    coolTemp: coolMean.get(bandName)
  };
};
// ========================================================================
// 1. LANDSAT 8/9 PROCESSING
// ========================================================================

print("Processing Landsat 8/9 data...");

// Load and process Landsat
var landsatCollection = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
  .merge(ee.ImageCollection("LANDSAT/LC09/C02/T1_L2"))
  .filterDate("2023-06-01", "2023-08-31")
  .filterBounds(losAngeles)
  .filter(ee.Filter.lt("CLOUD_COVER", 40))
  .map(maskLandsatClouds)
  .map(calculateLandsatLST);

// Create best-pixel composite
var landsatComposite = landsatCollection.qualityMosaic("LST");

// Calculate Landsat UHI
var landsatUHI = calculateUHI(landsatComposite, losAngeles, "LST");

print("Landsat Results:");
print("  Mean Temperature (°C):", landsatUHI.meanTemp);
print("  UHI Intensity (°C):", landsatUHI.uhiIntensity);

// ========================================================================
// 2. SENTINEL-2 PROCESSING
// ========================================================================

print("Processing Sentinel-2 data...");

// Load Sentinel-2
var sentinel2Collection = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterDate("2023-06-01", "2023-08-31")
  .filterBounds(losAngeles)
  .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 30))
  .map(processSentinel2);

// Create median composite
var sentinel2Composite = sentinel2Collection.median();

// Calculate Sentinel-2 UHI
var sentinelUHI = calculateUHI(sentinel2Composite, losAngeles, "Estimated_LST");

print("Sentinel-2 Results:");
print("  Mean Temperature (°C):", sentinelUHI.meanTemp);
print("  UHI Intensity (°C):", sentinelUHI.uhiIntensity);

// ========================================================================
// 3. SIMPLE LAND COVER ANALYSIS
// ========================================================================

// Classify areas based on NDVI only
var classifyByVegetation = function(image, ndviName) {
  var ndvi = image.select(ndviName);

  // Simple 3-class system based on vegetation
  var landCover = ee.Image.constant(0)
                    .where(ndvi.gt(0.6), 1)    // Dense vegetation
                    .where(ndvi.gt(0.3).and(ndvi.lte(0.6)), 2)  // Moderate vegetation
                    .where(ndvi.lte(0.3), 3)   // Low/no vegetation (urban)
                    .rename("LandCover");

  return landCover;
};

var landsatLandCover = classifyByVegetation(landsatComposite, "NDVI");
var sentinelLandCover = classifyByVegetation(sentinel2Composite, "NDVI");

// ========================================================================
// 4. VISUALIZATION
// ========================================================================

// Visualization parameters
var lstVis = {
  min: 25, max: 45,
  palette: ["#000080", "#0000FF", "#00FFFF", "#FFFF00", "#FF8000", "#FF0000"]
};

var uhiVis = {
  min: -5, max: 5,
  palette: ["#0000FF", "#80FFFF", "#FFFFFF", "#FF8080", "#FF0000"]
};

var ndviVis = {
  min: 0, max: 1,
  palette: ["#d7191c", "#fdae61", "#ffffbf", "#a6d96a", "#1a9641"]
};

var landCoverVis = {
  min: 0, max: 3,
  palette: ["#000000", "#228B22", "#90EE90", "#FF6B6B"]
};

// Center map on St. Louis
Map.centerObject(losAngeles, 10);

// Add Landsat layers
Map.addLayer(landsatComposite.select("LST").clip(losAngeles), lstVis, 
             "Landsat LST (°C)");
Map.addLayer(landsatUHI.uhiMap.clip(losAngeles), uhiVis, 
             "Landsat UHI Map", false);
Map.addLayer(landsatComposite.select("NDVI").clip(losAngeles), ndviVis, 
             "Landsat NDVI", false);

// Add Sentinel-2 layers
Map.addLayer(sentinel2Composite.select("Estimated_LST").clip(losAngeles), lstVis, 
             "Sentinel-2 LST (°C)", false);
Map.addLayer(sentinelUHI.uhiMap.clip(losAngeles), uhiVis, 
             "Sentinel-2 UHI Map", false);
Map.addLayer(sentinel2Composite.select("NDVI").clip(losAngeles), ndviVis, 
             "Sentinel-2 NDVI", false);

// Add land cover layers
Map.addLayer(landsatLandCover.clip(losAngeles), landCoverVis, 
             "Landsat Land Cover", false);
Map.addLayer(sentinelLandCover.clip(losAngeles), landCoverVis, 
             "Sentinel-2 Land Cover", false);

// ========================================================================
// 5. STATISTICS
// ========================================================================

// Land cover temperature analysis
var analyzeVegetationTemperature = function(lstImage, landCoverImage, sensorName) {
  print("\\n=== " + sensorName + " VEGETATION ANALYSIS ===");

  var lstBand = lstImage.bandNames().get(0);
  var classNames = ["Water/Other", "Dense Vegetation", "Moderate Vegetation", "Low Vegetation (Urban)"];

  for (var i = 0; i < 4; i++) {
    var mask = landCoverImage.eq(i);
    var temp = lstImage.updateMask(mask).reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: losAngeles,
      scale: 100,
      maxPixels: 1e9,
      bestEffort: true
    });

    var area = mask.multiply(ee.Image.pixelArea()).reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: losAngeles,
      scale: 100,
      maxPixels: 1e9,
      bestEffort: true
    });

    print("  " + classNames[i] + ":");
    print("    Temperature:", temp.get(lstBand));
    print("    Area (km²):", ee.Number(area.get("LandCover")).divide(1e6));
  }
};

// Run vegetation analysis for both sensors
analyzeVegetationTemperature(landsatComposite.select("LST"), landsatLandCover, "Landsat");
analyzeVegetationTemperature(sentinel2Composite.select("Estimated_LST"), sentinelLandCover, "Sentinel-2");

// Overall comparison
print("\\n=== SENSOR COMPARISON ===");
var landsatRange = landsatComposite.select("LST").reduceRegion({
  reducer: ee.Reducer.minMax(),
  geometry: losAngeles,
  scale: 100,
  maxPixels: 1e9
});

var sentinelRange = sentinel2Composite.select("Estimated_LST").reduceRegion({
  reducer: ee.Reducer.minMax(),
  geometry: losAngeles,
  scale: 100,
  maxPixels: 1e9
});

print("Landsat LST Range:", landsatRange);
print("Sentinel-2 LST Range:", sentinelRange);

// NDVI comparison
var landsatNDVI = landsatComposite.select("NDVI").reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: losAngeles,
  scale: 100,
  maxPixels: 1e9
});

var sentinelNDVI = sentinel2Composite.select("NDVI").reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: losAngeles,
  scale: 100,
  maxPixels: 1e9
});

print("Mean NDVI - Landsat:", landsatNDVI.get("NDVI"));
print("Mean NDVI - Sentinel-2:", sentinelNDVI.get("NDVI"));

//export function
var exportResults = false; // Set to true to enable

if (exportResults) {
  Export.image.toDrive({
    image: landsatComposite.select(["LST", "NDVI"]).clip(losAngeles),
    description: "Landsat_UHI_losAngeles_2023",
    scale: 100,
    region: losAngeles,
    maxPixels: 1e9
  });

  Export.image.toDrive({
    image: sentinel2Composite.select(["Estimated_LST", "NDVI"]).clip(losAngeles),
    description: "Sentinel2_UHI_losAngeles_2023",
    scale: 20,
    region: losAngeles,
    maxPixels: 1e9
  });
}
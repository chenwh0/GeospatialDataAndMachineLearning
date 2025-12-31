var start = "2023-04-01";
var end = "2023-10-31";


// Sentinel 2
var sentinel2_data = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(geometry)
  .filterDate(start, end)
  .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 40));
print("Sentinel-2 collection size:", sentinel2_data.size());

function mask_SCL_band(image) {
  var scl = image.select("SCL");
  var mask = scl.neq(3) // cloud shadow
  .and(scl.neq(8)) // medium clouds
  .and(scl.neq(9)) // high clouds
  .and(scl.neq(10)) // cirrus
  .and(scl.neq(11)) // snow/ice
  return image.updateMask(mask);
}


// create composite and export

var sentinel2composite = sentinel2_data.map(mask_SCL_band).median().clip(geometry);
var s2bands = ["B2", "B3", "B4", "B8"] // Blue, Green, Red, NIR
//Map.addLayer(s2composite, {bands: ["B4", "B3", "B2"], min: 0, max: 3000, gamma: 1.1}, "S2 true color");


var ndvi_sentinel2 = sentinel2composite.normalizedDifference(["B8", "B4"]).rename("NDVI");
var ndviVisualize = {min: -0.2, max: 0.8, palette: ["blue", "white", "green"]};
//Map.addLayer(ndvi_sentinel2, ndviVisualize, "S2 NDVI");

var savi_sentinel2 = sentinel2composite.expression("((NIR-RED) / (NIR+RED+L)) * (1+L)", {
 "NIR": sentinel2composite.select("B8"),
 "RED": sentinel2composite.select("B4"),
 "L": 0.5
}).rename("SAVI");
var saviVisualize = {min: -0.2, max: 0.6, palette: ["brown", "yellow", "green"]};
//Map.addLayer(savi_sentinel2, saviVisualize, "S2 SAVI");

var evi_sentinel2 = sentinel2composite.expression("2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))", {
 "NIR": sentinel2composite.select("B8"),
 "RED": sentinel2composite.select("B4"),
 "BLUE": sentinel2composite.select("B2")
}).rename("EVI");
var eviVisualize = {min: -0.2, max: 1.0, palette: ["blue", "white", "green"]};
//Map.addLayer(evi_sentinel2, eviVisualize, "S2 EVI");

var ndwi_sentinel2 = sentinel2composite.expression("(GREEN - NIR) / (GREEN + NIR)", {"GREEN": sentinel2composite.select("B3"), "NIR": sentinel2composite.select("B8")}).rename("NDWI");
var ndwiVisualize = {min: -0.3, max: 0.3, palette: ["red", "white", "blue"]};
//Map.addLayer(ndwi_sentinel2, ndwiVisualize, "S2 NDWI");


/*
// Export to google drive
Export.image.toDrive({
  image: s2composite.select(s2bands).toFloat(),
  description: "s2_FresnoCalifornia_2023",
  region: geometry,
  scale: 10,
  crs: "EPSG:4326", // For Fresno, California
  maxPixels: 1e13
});
*/



// Landsat 8 and 9
var landsat8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
  .filterBounds(geometry)
  .filterDate(start, end);
var landsat9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
  .filterBounds(geometry)
  .filterDate(start, end);

// Merge collections & apply cloud cover pre-filter
var landsat89 = landsat8.merge(landsat9).filter(ee.Filter.lt("CLOUD_COVER", 60));
print("Landsat8 and 9 collection size:", landsat89.size());

function quality_mask(image) {
  var quality = image.select("QA_PIXEL");
  
  // Define bit positions for quality flags
  var cloudShadowBitMask = 1 << 3;
  var snowBitMask = 1 << 4; 
  var cloudBitMask = 1 << 5;
  var cirrusBitMask = 1 << 7;
  
  // Combined mask
  var mask = quality.bitwiseAnd(cloudShadowBitMask).eq(0)
      .and(quality.bitwiseAnd(snowBitMask).eq(0))
      .and(quality.bitwiseAnd(cloudBitMask).eq(0))
      .and(quality.bitwiseAnd(cirrusBitMask).eq(0));
  
  var optical = image.select(["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"])
                      .multiply(0.0000275).add(-0.2);
  return image.addBands(optical, null, true).updateMask(mask);
}

var landsat89Masked = landsat89.map(quality_mask);
var seasonal_composite = landsat89Masked.median().clip(geometry); 
var bands = ["SR_B2", "SR_B3", "SR_B4", "SR_B5"]; // Blue, Green, Red, NIR
//Map.addLayer(seasonal_composite, {bands: ["SR_B4", "SR_B3", "SR_B2"], min: 0.02, max: 0.3, gamma: 1.1}, "Landsat true color");


var ndvi_landsat89 = seasonal_composite.normalizedDifference(['SR_B5','SR_B4']).rename('NDVI');
Map.addLayer(ndvi_landsat89, ndviVisualize, 'Landsat NDVI');

/*
Export.image.toDrive({
  image: seasonal_composite.select(bands).toFloat(),
  description: "Landsat8_and9_2023_Apr_Oct",
  region: geometry,
  scale: 30, // Native 30km resolution
  crs: "EPSG:4326", // For Fresno, California
  maxPixels: 1e13
});
*/
// Interested study area
var fire_region = ee.Geometry.Rectangle([-102,35.0, -99.5,36.0]);
Map.centerObject(fire_region, 9);

// 21 days before and after fire
var preStart = "2024-02-01";
var preEnd = "2024-02-25";
var postStart = "2024-03-16";
var postEnd = "2024-04-06";


function maskS2Clouds(img) {
  var scl = img.select('SCL');

  // bad classes
  var shadows = scl.eq(3);      // cloud shadow
  var clouds = scl.eq(8);       // cloud medium probability
  var cirrus = scl.eq(9);       // cirrus
  var snow = scl.eq(11);        // snow/ice

  var mask = shadows.or(clouds).or(cirrus).or(snow).not();

  return img.updateMask(mask);
}
// Sentinel-2 SR corrected
var s2 = ee.ImageCollection("COPERNICUS/S2_SR")
  .filterBounds(fire_region)
  .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 40))
  .select(['B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11','B12','SCL'])
  .map(maskS2Clouds);
/*
// Firgure 3
function calculateComposites(image) {
  var nbr = image.normalizedDifference(["B8", "B12"]).rename("NBR");
  var ndvi = image.normalizedDifference(["B8", "B4"]).rename("NDVI");
  var ndmi = image.normalizedDifference(["B8", "B11"]).rename("NDMI");
  var ndwi = image.normalizedDifference(["B3", "B8"]).rename('NDWI');

  var evi = image.expression(
    "2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))",
    {
      "NIR": image.select("B8").divide(10000),
      "RED": image.select("B4").divide(10000),
      "BLUE": image.select("B2").divide(10000)
    }
  ).rename("EVI");
  var red = image.select("B4").divide(10000);
  var nir = image.select("B8").divide(10000);
  var bai = ee.Image(1).divide(
    (ee.Image(0.1).subtract(red)).pow(2).add(
      (ee.Image(0.06).subtract(nir)).pow(2)
    )
  ).rename("BAI");
  
  var B04 = image.select('B4').divide(10000);
  var B06 = image.select('B6').divide(10000);
  var B07 = image.select('B7').divide(10000);
  var B8A = image.select('B8A').divide(10000);
  var B12 = image.select('B12').divide(10000);

  // Avoid divisions by zero and negative roots with small eps
  var eps = 1e-6;
  var term1 = B06.multiply(B07).multiply(B8A)
                .divide(B04.max(eps))
                .max(0).sqrt();
  var term2 = B12.subtract(B8A)
                .divide(B12.add(B8A).max(eps).sqrt())
                .add(1);

  var bais2 = ee.Image(1).subtract(term1).multiply(term2).rename('BAIS2');
  
  // MIRBI
  var swir1 = image.select("B11").divide(10000);
  var swir2 = image.select("B11").divide(10000);
  var mirbi = swir2.multiply(10)
    .subtract(swir1.multiply(0.98))
    .add(2)
    .rename("MIRBI");
  return image.addBands([nbr, ndvi, ndmi, ndwi, evi, bai, bais2, mirbi]);
}


var preFireComposite = s2.filterDate(preStart, preEnd)
  .map(calculateComposites)
  .median()
  .clip(fire_region);

var postFireComposite = s2.filterDate(postStart, postEnd)
  .map(calculateComposites)
  .median()
  .clip(fire_region);

var dNBR = preFireComposite.select("NBR").subtract(postFireComposite.select("NBR")).rename('dNBR');
var dNDVI = preFireComposite.select("NDVI").subtract(postFireComposite.select("NDVI")).rename("dNDVI");
var dNDMI = preFireComposite.select("NDMI").subtract(postFireComposite.select("NDMI")).rename("dNDMI");
var dEVI = preFireComposite.select("EVI").subtract(postFireComposite.select("EVI")).rename("dEVI");
var postBAI = postFireComposite.select("BAI").rename("BAI");
var postBAIS2 = postFireComposite.select("BAIS2").rename("BAIS2");
  
var compositeSeverity = dNBR.multiply(0.4)
  .add(dNDVI.multiply(0.2))
  .add(dNDMI.multiply(0.15))
  .add(dEVI.multiply(0.1))
  .add(postBAI.subtract(100).divide(1000).multiply(0.1))
  .add(postBAIS2.multiply(0.05))
  .rename("CompositeSeverity");
  
var severityNormalized = compositeSeverity.unitScale(-0.5, 2.0).rename("SeverityNormalized");

var recoveryPotential = preFireComposite.select("NDVI")
  .multiply(0.4)
  .add(preFireComposite.select("NDMI").multiply(0.3))
  .add(ee.Image(1).subtract(severityNormalized).multiply(0.3))
  .rename("RecoveryPotential");

/*
var recoveryParams = {min: 0, max: 1, palette: ['red', 'orange', 'yellow', 'lightgreen', 'darkgreen']};
Map.addLayer(recoveryPotential, recoveryParams, "Recovery Potential");


var fireConfidence = dNBR.gt(0.1)
  .and(dNDVI.gt(0.05))                               // NDVI shows vegetation loss
  .and(postBAI.gt(50))
  .rename("FireConfidence");

var confidenceParams = {min: 0, max: 1, palette: ['white', 'red']};
var recoveryParams = {min: 0, max: 1, palette: ['red', 'orange', 'yellow', 'lightgreen', 'darkgreen']};

Map.addLayer(fireConfidence, confidenceParams, "Fire Detection Confidence");


var severityClasses = severityNormalized
  .where(severityNormalized.lt(0.1), 0)                    // Unburned
  .where(severityNormalized.gte(0.1).and(severityNormalized.lt(0.3)), 1)   // Low
  .where(severityNormalized.gte(0.3).and(severityNormalized.lt(0.5)), 2)   // Moderate-Low
  .where(severityNormalized.gte(0.5).and(severityNormalized.lt(0.7)), 3)   // Moderate-High
  .where(severityNormalized.gte(0.7), 4)                   // High
  .rename('SeverityClasses')
  .toInt();

/*
// Visualization parameters
var compositeParams = {min: 0, max: 1, palette: ['green', 'yellow', 'orange', 'red', 'darkred']};
var classParams = {min: 0, max: 4, palette: ['green', 'yellow', 'orange', 'red', 'purple']};
var severityMask = severityClasses.updateMask(fireConfidence);

// Map layers
Map.addLayer(severityNormalized, compositeParams, 'Composite Severity (Normalized)');
Map.addLayer(severityMask, classParams, "Burn Area Only - Severity Classes");
*/

/*
// Save areaStats.csv
var pixelArea = ee.Image.pixelArea();

var mask0 = severityClasses.eq(0).multiply(fireConfidence);
var mask1 = severityClasses.eq(1).multiply(fireConfidence);
var mask2 = severityClasses.eq(2).multiply(fireConfidence);
var mask3 = severityClasses.eq(3).multiply(fireConfidence);
var mask4 = severityClasses.eq(4).multiply(fireConfidence);

function areaOf(mask, classID) {
  var area_m2 = mask.multiply(pixelArea).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: fire_region,
    scale: 20,
    maxPixels: 1e13
  }).get('SeverityClasses');   // band name doesn't matter, gets single value

  return ee.Feature(null, {
    SeverityClass: classID,
    Area_ha: ee.Number(area_m2).divide(10000)
  });
}

var fc = ee.FeatureCollection([
  areaOf(mask0, 0),
  areaOf(mask1, 1),
  areaOf(mask2, 2),
  areaOf(mask3, 3),
  areaOf(mask4, 4)
]);

print("Area Stats (Ha)", fc);

// Export CSV
Export.table.toDrive({
  collection: fc,
  description: 'LastName_FirstName_A501_AreaStats',
  fileFormat: 'CSV',
  fileNamePrefix: 'LastName_FirstName_A501_AreaStats'
});
*/



/*
var dNBR = preFireComposite.select("NBR").subtract(postFireComposite.select("NBR")).rename('dNBR');
var dNDVI = preFireComposite.select("NDVI").subtract(postFireComposite.select("NDVI")).rename("dNDVI");
var dNDMI = preFireComposite.select("NDMI").subtract(postFireComposite.select("NDMI")).rename("dNDMI");
var dEVI = preFireComposite.select("EVI").subtract(postFireComposite.select("EVI")).rename("dEVI");
var postBAI = postFireComposite.select("BAI").rename("BAI");
var postBAIS2 = postFireComposite.select("BAIS2").rename("BAIS2");
*/

/*
// Save to indexStats.csv
var indices = ["NBR", "NDVI", "NDMI", "NDWI", "EVI", "BAI", "BAIS2", "MIRBI"];
var indexStats = postFireComposite.select(['NBR','NDVI','NDMI','EVI','BAI','BAIS2'])
  .reduceRegion({
    reducer: ee.Reducer.mean()
             .combine(ee.Reducer.median(), '', true)
             .combine(ee.Reducer.min(), '', true)
             .combine(ee.Reducer.max(), '', true)
             .combine(ee.Reducer.stdDev(), '', true),
    geometry: fire_region,
    scale: 20,        // slightly coarser can help
    bestEffort: true, // GEE chooses a scale if too many pixels
    maxPixels: 1e13
  });

print(indexStats)

var indexStatsFC = ee.FeatureCollection([ee.Feature(null, indexStats)]);

Export.table.toDrive({
collection: indexStatsFC,
description: 'Chen_WenHsin_A501_IndexStats',
fileFormat: 'CSV',
fileNamePrefix: 'Chen_WenHsin_A501_IndexStats'
});
*/

/*
var compositeSeverity = dNBR.multiply(0.4)
  .add(dNDVI.multiply(0.2))
  .add(dNDMI.multiply(0.15))
  .add(dEVI.multiply(0.1))
  .add(postBAI.subtract(100).divide(1000).multiply(0.1))
  .add(postBAIS2.multiply(0.05))
  .rename("CompositeSeverity");
  
var severityNormalized = compositeSeverity.unitScale(-0.5, 2.0).rename("SeverityNormalized");

var recoveryPotential = preFireComposite.select("NDVI")
  .multiply(0.4)
  .add(preFireComposite.select("NDMI").multiply(0.3))
  .add(ee.Image(1).subtract(severityNormalized).multiply(0.3))
  .rename("RecoveryPotential");

var recoveryParams = {min: 0, max: 1, palette: ['red', 'orange', 'yellow', 'lightgreen', 'darkgreen']};
Map.addLayer(recoveryPotential, recoveryParams, "Recovery Potential");

/*

var fireConfidence = dNBR.gt(0.1)
  .rename("FireConfidence");

var confidenceParams = {min: 0, max: 1, palette: ['white', 'red']};
var recoveryParams = {min: 0, max: 1, palette: ['red', 'orange', 'yellow', 'lightgreen', 'darkgreen']};

Map.addLayer(fireConfidence, confidenceParams, "Fire Detection Confidence");


var severityClasses = severityNormalized
  .where(severityNormalized.lt(0.1), 0)                    // Unburned
  .where(severityNormalized.gte(0.1).and(severityNormalized.lt(0.3)), 1)   // Low
  .where(severityNormalized.gte(0.3).and(severityNormalized.lt(0.5)), 2)   // Moderate-Low
  .where(severityNormalized.gte(0.5).and(severityNormalized.lt(0.7)), 3)   // Moderate-High
  .where(severityNormalized.gte(0.7), 4)                   // High
  .rename('SeverityClasses')
  .toInt();


// Visualization parameters
var compositeParams = {min: 0, max: 1, palette: ['green', 'yellow', 'orange', 'red', 'darkred']};
var classParams = {min: 0, max: 4, palette: ['green', 'yellow', 'orange', 'red', 'purple']};
var severityMask = severityClasses.updateMask(fireConfidence);

// Map layers
Map.addLayer(severityNormalized, compositeParams, 'Composite Severity (Normalized)');
Map.addLayer(severityMask, classParams, "Burn Area Only - Severity Classes")

// Individual difference layers for comparison

/*
Map.addLayer(dNDVI, {min: -0.2, max: 0.8, palette: ['blue', 'white', 'red']}, 'dNDVI', false);
Map.addLayer(dNDMI, {min: -0.5, max: 1.0, palette: ['blue', 'white', 'red']}, 'dNDMI', false);
Map.addLayer(dEVI, {min: -0.1, max: 0.5, palette: ['blue', 'white', 'red']}, 'dEVI', false);

/*
// Figure 2
// MIRBI for Sentinel-2
function addMIRBI(img) {
  var swir2 = img.select("B12").divide(10000);
  var swir1 = img.select("B11").divide(10000);
  var mirbi = swir2.multiply(10)
    .subtract(swir1.multiply(0.98))
    .add(2)
    .rename("MIRBI");
  return img.addBands(mirbi);
}

var mirbiCollection = s2.map(addMIRBI);
var preMIRBI = mirbiCollection.filterDate(preStart, preEnd).median().select("MIRBI").clip(fire_region);
var postMIRBI = mirbiCollection.filterDate(postStart, postEnd).median().select("MIRBI").clip(fire_region);

var mirbiParams = {min: 0, max: 5, palette: ["white", "orange", "red", "black"]}
Map.addLayer(preMIRBI, mirbiParams, "Texas pre-fire MIRBI (Sentinel");
Map.addLayer(postMIRBI, mirbiParams, "Texas post-fire MIRBI (Sentinel2)");

/*
// BAIS2 for Sentinel-2
function addBAIS2(img) {
  var B04 = img.select("B4").divide(10000);
  var B06 = img.select("B6").divide(10000);
  var B07 = img.select("B7").divide(10000);
  var B8A = img.select("B8A").divide(10000);
  var B12 = img.select("B12").divide(10000);
  
  // Prevent division by 0 or negative roots
  var eps = 1e-6; // epsilon
  var term1 = B06.multiply(B07).multiply(B8A)
                .divide(B04.max(eps))
                .max(0).sqrt();
  var term2 = B12.subtract(B8A)
                .divide(B12.add(B8A).max(eps).sqrt())
                .add(1);
  var bais2 = ee.Image(1).subtract(term1).multiply(term2).rename("BAIS2");
  return img.addBands(bais2);
}

var bais2Collection = s2.map(addBAIS2);
var preBAIS2 = bais2Collection.filterDate(preStart, preEnd).median().select("BAIS2").clip(fire_region);
var postBAIS2 = bais2Collection.filterDate(postStart, postEnd).median().select("BAIS2").clip(fire_region);

var bais2Params = {min: -1, max: 1, palette: ["white", "purple", "black"]};
Map.addLayer(preBAIS2, bais2Params, "Texas pre-fire BAIS2 (Sentinel2)");
Map.addLayer(postBAIS2, bais2Params, "Texas post-fire BAIS2 (Sentinel2)");
*/
/*
// BAI for Sentinel2
function addBAI(img) {
  var red = img.select("B4").divide(10000);
  var nir = img.select("B8").divide(10000);
  var bai = ee.Image(1).divide(
    (red.subtract(0.1)).pow(2).add((nir.subtract(0.06)).pow(2))
  ).rename("BAI");
  return img.addBands(bai);
}

var baiCollection = s2.map(addBAI);
var preBAI = baiCollection.filterDate(preStart, preEnd).median().select("BAI").clip(fire_region);
var postBAI = baiCollection.filterDate(postStart, postEnd).median().select("BAI").clip(fire_region);
var baiParams = {min: 0, max: 800, palette: ["white", "black"]};
//Map.addLayer(preBAI, baiParams, "Texas pre-fire BAI (Sentinel-2)");
//Map.addLayer(postBAI, baiParams, "Texas post-fire BAI (Sentinel-2)");
*/




// Figure 1
// 4) Function to calculate NBR for Sentinel-2
// NBR = (NIR - SWIR2) / (NIR + SWIR2) => B8 and B12 for Sentinel-2
function addNBR(image) {
  var nbr = image.normalizedDifference(["B8", "B12"]).rename("NBR");
  return image.addBands(nbr);
}

// 5) Apply to image collection
var nbrCollection = s2.map(addNBR);

// 6) Build pre-fire and post-fire median NBR composites over California
var preFireNBR = nbrCollection
  .filterDate(preStart, preEnd)
  .median()
  .select("NBR")
  .rename("preNBR")
  .clip(fire_region);

var postFireNBR = nbrCollection
  .filterDate(postStart, postEnd)
  .median()
  .select("NBR")
  .rename("postNBR")
  .clip(fire_region);

// 7) Optional: compute dNBR and a simple severity view
var dNBR = preFireNBR.subtract(postFireNBR).rename("dNBR");
var severity = dNBR.multiply(1000)
  .where(dNBR.multiply(1000).lt(100), 0)                            // unburned/very low
  .where(dNBR.multiply(1000).gte(100).and(dNBR.multiply(1000).lt(270)), 1)   // low
  .where(dNBR.multiply(1000).gte(270).and(dNBR.multiply(1000).lt(440)), 2)   // moderate-low
  .where(dNBR.multiply(1000).gte(440).and(dNBR.multiply(1000).lt(660)), 3)   // moderate-high
  .where(dNBR.multiply(1000).gte(660), 4)                           // high
  .clip(fire_region);

// 8) Visualization parameters
var nbrParams = {
  min: -1.0,
  max: 1.0,
  palette: ["red", "orange", "yellow", "white", "cyan", "blue"]
};

var dnbrParams = {min: -0.2, max: 0.8, palette: ["blue", "white", "red"]};
var sevParams  = {min: 0, max: 4, palette: ["green", "yellow", "orange", "red", "purple"]};

// 9) Map layers
Map.addLayer(preFireNBR, nbrParams, "Texas Pre-Fire NBR");
Map.addLayer(postFireNBR, nbrParams, "Texas Post-Fire NBR");
Map.addLayer(dNBR, dnbrParams, "Texas dNBR");
//Map.addLayer(severity, sevParams, "Texas Burn Severity (simple classes)");
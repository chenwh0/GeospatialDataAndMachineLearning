var startDate = "2023-06-01";
var endDate = "2023-09-30";

var region = ee.Geometry.Rectangle([120.18,23.57, 121.08,23.77]);


// Cloud-masking function
function maskL8(image) {
  var qa = image.select('QA_PIXEL');
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 4);
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
               .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}

// Data
var landsat8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
              .filterBounds(region)
              .filterDate(startDate, endDate)
              .filter(ee.Filter.lt("CLOUD_COVER", 10))
              .map(maskL8);  // Apply cloud-masking
              
var landsat8_composite = landsat8.median().clip(region);

function scaleLandsat(img) {
  var bands = ["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"];
  var scaled = img.select(bands).multiply(0.0000275).add(-0.2);
  return img.addBands(scaled, bands, true)
}


var landsat_image = scaleLandsat(landsat8_composite);
var ndvi = landsat_image.normalizedDifference(["SR_B5", "SR_B4"]).rename("NDVI");    
var ndwi = landsat_image.normalizedDifference(["SR_B3", "SR_B5"]).rename("NDWI");
landsat_image = landsat_image.addBands([ndvi, ndwi]);

var predictors = landsat_image.select(["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7", "NDVI", "NDWI"]);



// Prepare WorldCover labels
var wc = ee.Image('ESA/WorldCover/v200/2021');
var wcCodes = [10,20,30,40,50,60,70,80,90,95,100];
var wcNames = [
  'Tree cover','Shrubland','Grassland','Cropland','Built-up',
  'Bare / sparse','Snow & ice','Water','Herbaceous wetland',
  'Mangroves','Moss & lichen'
];
var wcPalette = [
  '006400','ffbb22','ffff4c','f096ff','fa0000',
  'b4b4b4','f0f0f0','0064c8','0096a0','00cf75','fae6a0'
];

var wcCodesList = ee.List(wcCodes);
var toIdx = ee.List.sequence(0, wcCodes.length - 1);
var ones = ee.List.repeat(1, wcCodesList.length());
var keptMask = wc.select('Map').remap(wcCodesList, ones).neq(0);

var labels = wc.select('Map')
  .updateMask(keptMask)
  .remap(wcCodesList, toIdx)
  .rename('landcover')
  .toInt();
  
  
// Prepare training & test data
var sample = predictors.addBands(labels).stratifiedSample({
  numPoints: 3000,
  classBand: "landcover",
  region: region,
  scale: 10,
  geometries: true,
  seed: 42
});

var withRand = sample.randomColumn("rand", 42); // Add random number between 0 & 1 column to every sample
var train = withRand.filter(ee.Filter.lt("rand", 0.7)); // Values > 0.7 will become training samples
var test = withRand.filter(ee.Filter.gte("rand", 0.7)); // Rest will be test samples

// 3 Classifiers
// Random forest
var random_forest = ee.Classifier.smileRandomForest({
  numberOfTrees: 300,
  seed: 42
}).train({
  features: train,
  classProperty: "landcover",
  inputProperties: predictors.bandNames()
});

// Support Vector Machine (SVM)
var support_vector_machine = ee.Classifier.libsvm({
  kernelType: "RBF",
  gamma: 0.5,
  cost: 10
}).train({
  features: train,
  classProperty: "landcover",
  inputProperties: predictors.bandNames()
});

// Classification And Regression Trees (CART)
var cart = ee.Classifier.smileCart({
  maxNodes: 100
}).train({
  features: train,
  classProperty: "landcover",
  inputProperties: predictors.bandNames()
});

// Execute classifiers
var random_forest_predictions = predictors.classify(random_forest).rename("RF");
var support_vector_machine_predictions = predictors.classify(support_vector_machine).rename("SVM");
var cart_predictions = predictors.classify(cart).rename("CART");

// Evaluation
function evaluateClassifier(classifier, classifer_name) {
  var test_predictions = test.classify(classifier);
  var confusion_matrix = test_predictions.errorMatrix("landcover", "classification");
  
  print("RESULTS for" + classifer_name);
  print("Confusion matrix:", confusion_matrix);
  print("Overall accuracy:", confusion_matrix.accuracy());
  print("Kappa:", confusion_matrix.kappa());
  print("Producers accuracy:", confusion_matrix.producersAccuracy());
  print("Users accuracy:", confusion_matrix.consumersAccuracy());
  
  return {
    name: classifer_name,
    accuracy: confusion_matrix.accuracy(),
    kappa: confusion_matrix.kappa(),
    confusion_matrix: confusion_matrix
  };
}

var random_forest_evaluation = evaluateClassifier(random_forest, "Random Forest");
var support_vector_machine_evaluation = evaluateClassifier(support_vector_machine, "Support Vector Machine");
var cart_evaluation = evaluateClassifier(cart, "CART");

// Calculate area statisics
function calculateAreaStats(predictions, classifier_name) {
  var area_image = ee.Image.pixelArea().addBands(predictions);
  var area_stats = area_image.reduceRegion({
    reducer: ee.Reducer.sum().group({
      groupField: 1,
      groupName: "class"
    }),
    geometry: region,
    scale: 30,
    maxPixels: 1e10
  });
  print("AREA STATS for", classifier_name);
  print("Area by class", area_stats);
  return area_stats;
}

var random_forest_stats = calculateAreaStats(random_forest_predictions, "Random Forest");
var support_vector_machine_stats = calculateAreaStats(support_vector_machine_predictions, "Support Vector Machine");
var cart_stats = calculateAreaStats(cart_predictions, "CART");

// Visualize
function buildLegend(title, labels, palette, showIndex) {
  var panel = ui.Panel({
    style: {position: "bottom-right",
            padding: "8px",
            backgroundColor: "white"}
  });
  panel.add(ui.Label({value: title, style: {fontWeight: "bold"}}));
  for (var i = 0; i < labels.length; i++) {
    var colorBox = ui.Label('', {
      backgroundColor: palette[i],
      padding: '8px',
      margin: '0 8px 4px 0'
    });
    var text = showIndex ? (i + ' — ' + labels[i]) : labels[i];
    var row = ui.Panel({
      widgets: [colorBox, ui.Label(text)],
      layout: ui.Panel.Layout.Flow('horizontal')
    });
    panel.add(row);
  }
  return panel;
}

var true_landsat_colors = {bands: ["SR_B4", "SR_B3", "SR_B2"], min: 0.05, max: 0.35, gamma: 1.2};
var wc_region = wc.select("Map").clip(region);
var wc_visualization = {min: 10, max: 100, palette: wcPalette};
var predicted_palette = wcPalette.slice();
var predicted_visualization = {min: 0, max: 10, palette: predicted_palette};
var legend = buildLegend("Classification legend", wcNames, predicted_palette, true);

// Map
Map.addLayer(region, {color: "pink"}, "study region", false);
Map.centerObject(region, 10);
Map.addLayer(landsat_image, true_landsat_colors, "True Landsat8");
Map.addLayer(wc_region, wc_visualization, "WorldCover (reference)");
Map.addLayer(random_forest_predictions, predicted_visualization, "Random Forest predictions");
Map.addLayer(support_vector_machine_predictions, predicted_visualization, "Support Vector Machine predictions");
Map.addLayer(cart_predictions, predicted_visualization, "CART predictions");
Map.add(legend);

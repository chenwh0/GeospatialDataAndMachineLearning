var coordinates = ee.Geometry.Point([-119.417931, 36.778259]); // Fresno, California
var stateName = "California"
var state = ee.FeatureCollection("TIGER/2018/States")
  .filter(ee.Filter.eq("NAME", stateName));
var stateGeometry = state.geometry();
// Time period for 2021-2022 drought analysis
var startDate = "2021-01-01";
var endDate = "2023-01-01";

// Get GRIDMET drought indices
var drought = ee.ImageCollection("GRIDMET/DROUGHT")
  .filterDate(startDate, endDate)
  .filterBounds(state);
  
// Extract time series for SPI, SPEI, EDDI, and PDSI
var extractDroughtIndices = function(image) {
  var date = image.date().format("YYYY-MM-dd");
  var values = image.select([
    "spi90d", // agricultural drought
    "spei90d", // comprehensive drought
    "eddi90d", // atmospheric demand
    "pdsi" // long-term drought
    ]).reduceRegion({
      reducer: ee.Reducer.first(),
      geometry: coordinates, 
      scale: 4638
    });
    return ee.Feature(null, {
      "date": date,
      "SPI_90d": values.get("spi90d"),
      "SPEI_90d": values.get("spei90d"),
      "EDDI_90d": values.get("eddi90d"),
      "PDSI": values.get("pdsi")
    });
};

// Time series feature collection
var droughtTimeSeries = drought.map(extractDroughtIndices);
var SPISeries = droughtTimeSeries.map(function(f) {
  return ee.Feature(null, {
    "date": f.get("date"),
    "value": f.get("SPI_90d"),
    "series": "SPI-90d"
  });
});
var SPEISeries = droughtTimeSeries.map(function(f) {
  return ee.Feature(null, {
    "date": f.get("date"),
    "value": f.get("SPEI_90d"),
    "series": "SPEI-90d"
  });
});
var EDDISeries = droughtTimeSeries.map(function(f) {
  return ee.Feature(null, {
    "date": f.get("date"),
    "value": f.get("EDDI_90d"),
    "series": "EDDI-90d"
  });
});
var PDSISeries = droughtTimeSeries.map(function(f) {
  return ee.Feature(null, {
    "date": f.get("date"),
    "value": f.get("PDSI"),
    "series": "PDSI"
  });
});

// Combine all series
var allSeries = SPISeries.merge(SPEISeries).merge(EDDISeries).merge(PDSISeries);

// Create comprehensive drought index chart
var droughtChart = ui.Chart.feature.groups(allSeries, "date", "value", "series")
  .setChartType("LineChart")
  .setOptions({
    title: "Multi-Index Drought Analysis — Fresno, California — 2021-2022",
    hAxis: {title: "Date"},
    vAxis: {title: "Index Value"},
    legend: {position: "top"},
    series: {
      0: {color: "#1f77b4", lineWidth: 2}, // SPI
      1: {color: "#ff7f0e", lineWidth: 2}, // SPEI  
      2: {color: "#2ca02c", lineWidth: 2}, // EDDI
      3: {color: "#d62728", lineWidth: 2}  // PDSI
    },
    interpolateNulls: true
  });

print(droughtChart);

// Add drought severity reference lines
print("Drought Severity Classifications:");
print("SPI/SPEI/EDDI: -0.8 (Mild), -1.3 (Moderate), -1.6 (Severe), -2.0 (Extreme)");
print("PDSI: -1.0 (Mild), -2.0 (Moderate), -3.0 (Severe), -4.0 (Extreme)");







// Spatial Drought mapping ()
var month = drought.filterDate("2021-06-01", "2021-07-01");
var monthDrought = ee.Algorithms.If(
  month.size().gt(0),
  month.sort("system:time_start", false).first(),
  drought.sort("system:time_start", false).first()
  );
monthDrought = ee.Image(monthDrought);

// indices for drought analysis
var spi90 = monthDrought.select("spi90d").clip(stateGeometry);
var spei90 = monthDrought.select("spei90d").clip(stateGeometry);
var pdsi = monthDrought.select("pdsi").clip(stateGeometry);

// Visualize
var droughtColors = [
  "#730000","#e60000","#ffaa00","#fcd37f","#ffff00",
  "#ffffff","#aaff55","#00ffff","#00aaff","#0000ff","#0000aa"
];
var standardVisualization = {min: -3, max: 3, palette: droughtColors};
var pdsiVisualization = {min: -6, max: 6, palette: droughtColors};

Map.centerObject(state, 6);
Map.addLayer(spi90, standardVisualization, "SPI-90d (June 2021)");
Map.addLayer(spei90, standardVisualization, "SPEI-90d (June 2021)", false);
Map.addLayer(pdsi, pdsiVisualization, "PDSI (June 2021)", false);
Map.addLayer(coordinates, {color: "black"}, "Fresno, California");


var agricultureVisualization = {min: -2.5, max: 2.5, palette: droughtColors}; // Agricultural focus
var springPlanting = drought.filterDate('2021-04-01', '2021-05-31');
var summerGrowth = drought.filterDate('2021-06-01', '2021-08-31'); 
var fallHarvest = drought.filterDate('2021-09-01', '2021-10-31');

// Function to calculate agricultural drought metrics
var calculateAgDrought = function(collection, period) {
  var agIndices = collection.select(['spi30d', 'spi90d', 'spei30d', 'spei90d']);
  
  var seasonalMean = agIndices.mean().set('period', period);
  
  var droughtFreq = agIndices.select('spi30d')
    .map(function(img) {
      return img.lt(-0.8);
    }).mean().multiply(100)
    .rename('drought_frequency_pct');
    
  return seasonalMean.addBands(droughtFreq);
};

// Calculate agricultural drought for each period
var springAg = calculateAgDrought(springPlanting, 'Spring_2022');
var summerAg = calculateAgDrought(summerGrowth, 'Summer_2022');
var fallAg = calculateAgDrought(fallHarvest, 'Fall_2022');

// ===============================
// Enhanced Map Visualizations
// ===============================

Map.centerObject(state, 6);

// 1) Seasonal drought progression maps using agVis parameters
Map.addLayer(springAg.select('spi30d').clip(stateGeometry), agricultureVisualization, 
            'Spring Agricultural Drought (SPI-30d)');
Map.addLayer(summerAg.select('spi30d').clip(stateGeometry), agricultureVisualization,
            'Summer Agricultural Drought (SPI-30d)');
Map.addLayer(fallAg.select('spi30d').clip(stateGeometry), agricultureVisualization,
            'Fall Agricultural Drought (SPI-30d)', false);

// 2) Multi-timescale comparison using stdVis parameters
Map.addLayer(summerAg.select('spi30d').clip(stateGeometry), standardVisualization,
            'Summer SPI-30d (Shorter-term)', false);
Map.addLayer(summerAg.select('spi90d').clip(stateGeometry), standardVisualization,
            'Summer SPI-90d (Longer-term)', false);


Map.addLayer(summerAg.select('spei30d').clip(stateGeometry), agricultureVisualization,
            'Summer SPEI-30d (Evapotranspiration)', false);

// 3) Drought frequency visualization
var freqVis = {min: 0, max: 100, palette: ['#00ff00','#ffff00','#ff8000','#ff0000','#800000']};
Map.addLayer(summerAg.select('drought_frequency_pct').clip(stateGeometry), freqVis,
            'Summer Drought Frequency (%)', false);
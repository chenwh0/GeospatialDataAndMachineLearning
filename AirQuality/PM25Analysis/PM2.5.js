// 2021-2022 daily PM2.5 time series chart from MERRA-2 data 
var interested_region = ee.Geometry.Point([-119.02, 35.37])
var start = "2021-01-01";
var end = "2022-12-31";
// PM2.5 reconstruction formula = Dust surface + Sea salt + Organic carbon + Black Carbon + Sulfate factor*SO4
var hourly_concentration = ee.ImageCollection("NASA/GSFC/MERRA/aer/2")
.select(["DUSMASS25", "SSSMASS25", "OCSMASS", "BCSMASS", "SO4SMASS"])
.filterDate(start, end)
.filterBounds(interested_region);

// Build in PM2.5 units kg/m^3
var toPM25 = function(img) {
  var so4_conv = img.select("SO4SMASS").multiply(1.375); // Sulfate factor
  var pm25 = img.select(["DUSMASS25", "SSSMASS25", "OCSMASS", "BCSMASS"])
  .addBands(so4_conv.rename("SO4_conv"))
  .reduce(ee.Reducer.sum())
  .rename("PM25"); 
  return pm25.copyProperties(img, ["system:time_start","system:time_end"]);
};

var hourlyPM = hourly_concentration.map(toPM25).filter(ee.Filter.notNull(["system:time_start"]));

var nDays = ee.Date(end).difference(ee.Date(start), "day");
var days = ee.List.sequence(0, nDays.subtract(1));

// Calculate daily kg/m^3
var dailyPM = ee.ImageCollection.fromImages(days.map(function(d){
  var d0 = ee.Date(start).advance(d, "day");
  var d1 = d0.advance(1, "day");
  var subset = hourlyPM.filterDate(d0, d1);
  return ee.Algorithms.If(
    subset.size().gt(0),
    subset.mean().set("system:time_start", d0.millis()),
    null
    );
  })).filter(ee.Filter.notNull(["system:time_start"]));

// Convert daily averages to µg/m^3 (microgram/m^3)
var dailyPM_micrograms = dailyPM.map(function(img) {
  var time_start = ee.Date(img.get("system:time_start"));
  var time_end = ee.Date(img.get("system:time_end"));
    return img.multiply(1e9).rename("PM25_microgramm3")
    .copyProperties(img,  ["system:time_start","system:time_end"]);
});

// Chart
var chart = ui.Chart.image.series({
  imageCollection: dailyPM_micrograms,
  region: interested_region,
  reducer: ee.Reducer.mean(),
  scale: 50000
}).setOptions({
  title: "Bakersfield, CA - Daily Surface PM2.5 (MERRA-2, 2021-2022)",
  vAxis: {title: "PM2.5 (µg/m^3)"},
  hAxis: {title: "Date"},
  interpolateNulls: true,
  lineWidth: 2,
  pointSize: 1,
});

print(chart);








// 2021 monthly aggregated PM2.5 chart from GHAP data
var start = "2021-01-01";
var end = "2021-12-31";

var ghap = ee.ImageCollection("projects/sat-io/open-datasets/GHAP/GHAP_D1K_PM25")
  .filterDate(start, end)
  .filterBounds(interested_region);

// Check availability
print("Total daily images:", ghap.size());

// Scale factor: 0.1 -> convert to µg/m³
var applyScale = function(img) {
  return img.multiply(0.1).copyProperties(img, ["system:time_start"]);
};
var daily = ghap.map(applyScale);

// Monthly means (12 images, tag month and canonical time)
var monthlyAverages = ee.ImageCollection.fromImages(
  ee.List.sequence(1, 12).map(function(m) {
    var i = daily.filter(ee.Filter.calendarRange(m, m, "month")).mean();
    return i.set("system:time_start", ee.Date.fromYMD(2021, m, 1));
  })
);

// Monthly time series chart
var chart = ui.Chart.image.series({
  imageCollection: monthlyAverages,
  region: interested_region,
  reducer: ee.Reducer.mean(),
  scale: 1000
}).setOptions({
  title: "Bakersfield, California - Monthly PM2.5 (GHAP, 2021)",
  vAxis: {title: "PM2.5 (µg/m³)"},
  hAxis: {title: "Date"},
  lineWidth: 3,
  pointSize: 4
});
print(chart);

// Optional: daily chart (can be heavy)
var dailyChart = ui.Chart.image.series({
  imageCollection: daily,
  region: roi,
  reducer: ee.Reducer.mean(),
  scale: 1000
}).setOptions({
  title: "Bakersfield, California - Daily PM2.5 (GHAP, 2020)",
  vAxis: {title: "PM2.5 (µg/m³)"},
  hAxis: {title: "Date"},
  lineWidth: 1
});
print(dailyChart);






// 2021 seasonal PM2.5 average comparison chart from GHAP data
// GHAP daily PM2.5 (1 km) -> apply scale factor 0.1 and rename to PM25
var daily = ee.ImageCollection("projects/sat-io/open-datasets/GHAP/GHAP_D1K_PM25")
  .filterDate(start, end)
  .filterBounds(interested_region)
  .map(function(img){
    return img.multiply(0.1)
              .rename("PM25")
              .copyProperties(img, ["system:time_start"]);
  });

// Monthly means (12 images, tag month and canonical time)
var monthlyAverages = ee.ImageCollection.fromImages(
  ee.List.sequence(1, 12).map(function(m){
    var mImg = daily.filter(ee.Filter.calendarRange(m, m, "month")).mean();
    return mImg.set({
      "month": m,
      "system:time_start": ee.Date.fromYMD(2020, m, 1).millis()
    });
  })
);

var addSeason = function(img){
  var m = ee.Number(img.get("month"));
  var season = ee.String(ee.Algorithms.If(
    m.gte(3).and(m.lte(5)), "Spring",
    ee.Algorithms.If(
      m.gte(6).and(m.lte(8)), "Summer",
      ee.Algorithms.If(
        m.gte(9).and(m.lte(11)), "Fall", "Winter"
      )
    )
  ));
  return img.set("season", season);
};
var monthlySeasoned = monthlyAverages.map(addSeason);

// Compute seasonal averages at ROI (returns a compact FeatureCollection)
var seasonNames = ee.List(["Spring", "Summer", "Fall", "Winter"]);
var seasonalFC = ee.FeatureCollection(seasonNames.map(function(s){
  var meanImg = monthlySeasoned.filter(ee.Filter.eq("season", s)).mean();
  var seasonAverage = meanImg.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: interested_region,
    scale: 1000
  }).get("PM25");
  return ee.Feature(null, {season: s, pm25: seasonAverage});
}));

print("Seasonal PM2.5 Averages (GHAP 2021, µg/m³):", seasonalFC);

// Column chart with clean axes (y starts at 0, x labels hidden)
var chart = ui.Chart.feature.byFeature({
  features: seasonalFC,
  xProperty: "season",
  yProperties: ["pm25"]
})
.setChartType("ColumnChart")
.setOptions({
  title: "Seasonal PM2.5 at Bakersfield, California (GHAP 2020)",
  hAxis: {title: "Season", ticks: []},
  vAxis: {title: "PM2.5 (µg/m³)", viewWindow: { min: 0 }    // start y-axis at 0
  },
  legend: { position: "none" },
  colors: ["#66c2a5"],
  bar: { groupWidth: "80%" }
});
print(chart);
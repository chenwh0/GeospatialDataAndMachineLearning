// Global vars
var interested_region = ee.Geometry.Point([-118.2437, 34.0522]); // Los Angeles, California
var start = "2024-01-01";
var end = "2024-12-31";
var no2 = ee.ImageCollection("COPERNICUS/S5P/OFFL/L3_NO2")
            .filterDate(start, end)
            .filterBounds(interested_region)
            .select("tropospheric_NO2_column_number_density"); // mol/m²




// Monthly Time Series (Sentinel-5P OFFL)
var months = ee.List.sequence(1,12);
var monthly = ee.ImageCollection.fromImages(months.map(function(month) {
        var monthAvg = no2.filter(ee.Filter.calendarRange(month, month, "month")).mean();
        return monthAvg.set({
            "system:time_start": ee.Date.fromYMD(2024, month, 1).millis(),
            "month": month
        });
    }
));

var chart = ui.Chart.image.series({
    imageCollection: monthly, 
    region: interested_region,
    reducer: ee.Reducer.mean(),
    scale: 11100 // S5P L3 optimal scale (~0.1 deg binning)
}).setOptions({
    title: "Los Angeles, California Monthly Tropospheric NO₂ averages (S5P OFFL, 2024)",
    hAxis: {title: "Month"},
    vAxis: {
        title: "NO₂ Column Density (mol/m²)",
        viewWindow: {min: 0},
        format: "scientific" // Better for scientific notation
    },
    lineWidth: 3,
    pointSize: 4,
    colors: ["#d62728"]
});
print(chart);
            
            
            
// Seasonal Pattern Analysis (Sentinel-5P OFFL) 

// Classify meteorological seasons
var addSeason = function(img) {
    // Extract month from system:time_start
    var date = ee.Date(img.get("system:time_start"));
    var month = date.get("month");

    var season =
        ee.Algorithms.If(month.gte(3).and(month.lte(5)), "Spring",
        ee.Algorithms.If(month.gte(6).and(month.lte(8)), "Summer",
        ee.Algorithms.If(month.gte(9).and(month.lte(11)), "Fall",
        "Winter")
        )
    );
    return img.set("season", season);
};

var dailySeasoned = no2.map(addSeason); // Apply seasonal classification to daily data
var seasonNames = ee.List(["Spring", "Summer", "Fall", "Winter"]);
var seasonalFeatures = ee.FeatureCollection(seasonNames.map(function(seasonName) {
    var seasonAvg = dailySeasoned
        .filter(ee.Filter.eq("season", seasonName))
        .mean();
    var seasonalValue = seasonAvg.reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: interested_region,
        scale: 11100, // S5P optimal scale
        maxPixels: 1e9
    }).get("tropospheric_NO2_column_number_density");
    return ee.Feature(null, {
        season: seasonName,
        no2: seasonalValue
    });

}));
print("Seasonal NO₂ averages (Los Angeles, California, 2024)", seasonalFeatures);

// Styled column with scientific format
var seasonalChart = ui.Chart.feature.byFeature({
    features: seasonalFeatures,
    xProperty: "season",
    yProperties: ["no2"]
})
.setChartType("ColumnChart")
.setOptions({
    title: "Los Angeles, California Seasonal Tropospheric NO₂ averages (S5P OFFL, 2024)",
    hAxis: {
        title: "Season",
        ticks: [] // clean x-axis, no tick marks
    },
    vAxis: {
        title: "NO₂ column density (mol/m²)",
        viewWindow: {min: 0},
        format: "scientific"
    },
    legend: {position: "none"},
    colors: ["#d62728"], // red for NO₂
    bar: {groupWidth: "80%"},
    backgroundColor: "white"
});
print(seasonalChart);
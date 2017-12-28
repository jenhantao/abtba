$(document).ready(function(){
 
    // Function to get the max value in an Array
    Array.max = function(array){
        return Math.max.apply(Math,array);
    };
 
    // Get all data values from our table cells making sure to ignore the first column of text
    // Use the parseInt function to convert the text string to a number
 
    var counts= $('.heat-map tbody td').not('.stats-title').map(function() {
        return parseInt($(this).text());
    }).get();
 
    // run max value function and store in variable
    var max = Array.max(counts);
 
    n = 100; // Declare the number of groups
 
    // Define the ending colour, which is white
    xr = 255; // Red value
    xg = 255; // Green value
    xb = 255; // Blue value
 
    // Define the starting colour #f32075
    yr = 243; // Red value
    yg = 32; // Green value
    yb = 117; // Blue value
 
    // Loop through each data point and calculate its % value
    $('.heat-map tbody td').not('.stats-title').each(function(){
        var val = parseInt($(this).text());
        var pos = parseInt((Math.round((val/max)*100)).toFixed(0));
        red = parseInt((xr + (( pos * (yr - xr)) / (n-1))).toFixed(0));
        green = parseInt((xg + (( pos * (yg - xg)) / (n-1))).toFixed(0));
        blue = parseInt((xb + (( pos * (yb - xb)) / (n-1))).toFixed(0));
        clr = 'rgb('+red+','+green+','+blue+')';
	$(this).css({backgroundColor:clr});
    });
});

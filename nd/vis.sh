#dim=$(((dim - 100)/$#))
var=$1
if [ "$var" -eq "$var" ] 2>/dev/null; then filter=$1;dm=$(($#-1));shift; else filter=0;dm=$#; fi

dim1=`xdpyinfo|grep dimensions| awk -F'[:xp]' '{print $3}'` #height
dim2=`xdpyinfo|grep dimensions| awk -F'[:x]' '{print $2}'|awk '{print $1}'` #width
dim1=$((dim1-100))
dim2=$(((dim2-100)/$#))
dim=$dim1
if [[ $dim2 -lt $dim1 ]]
  then  dim=$dim2
fi

tmp="tmp.html"
echo "" >$tmp
echo "
<!DOCTYPE html>
<meta charset=\"utf-8\">
<style>

.node {
  cursor: pointer;
  stroke: #000;
  stroke-width: 0.1px;  
}

.node:hover {
  stroke: #000 !important;
  stroke-width: 8px !important;
}

.node--leaf {
  fill: white;
}

.label {
  font: 16px \"Helvetica Neue\", Helvetica, Arial, sans-serif;
  text-anchor: middle;
  text-shadow: 0 1px 0 #fff, 1px 0 0 #fff, -1px 0 0 #fff, 0 -1px 0 #fff;
}

.label,
.node--root {
  pointer-events: none;
}

</style>
<body>
    <div class=\"container\" style=\"width: 100%;\">" >>$tmp
    
i=1;
for arg in "$@"
do
    echo "
    <div id=\""dv_$i"\" style=\"width: "$dim"px; float: left;\">
        <center>"$arg"</center>
    </div>" >>$tmp;
    i=$((i+1));
done

echo "
</div>  
</body>

<script src=\"http://d3js.org/d3.v3.min.js\"></script>
<script>

var margin = 20,
    diameter = "$dim";

var color = d3.scale.linear()
    .domain([0, 1, 2, 3, 4])
    .range([\"rgb(0,0,255)\", \"rgb(0,255,255)\", \"rgb(0,255,0)\", \"rgb(255,255,0)\", \"rgb(255,0,0)\"])
    .interpolate(d3.interpolate);

var pack = d3.layout.pack()
    .padding(2)
    .size([diameter - margin, diameter - margin])
    .value(function(d) { return d.size; })" >>$tmp
    
i=1;
for arg in "$@"
do
    echo "
var svg"$i" = d3.select(\"#"dv_$i"\").append(\"svg\")
    .attr(\"width\", diameter)
    .attr(\"height\", diameter)
    .append(\"g\")
    .attr(\"transform\", \"translate(\" + diameter / 2 + \",\" + diameter / 2 + \")\");" >>$tmp
    i=$((i+1));
done



i=1;
for arg in "$@"
do
    echo "
d3.json(\""$arg"\", function(error, root) {
  if (error) throw error;

  var focus = root,
      nodes = pack.nodes(root),
      view;

  var circle = svg"$i".selectAll(\"circle\")
      .data(nodes)      
    .enter().append(\"circle\")      
      .attr(\"class\", function(d) { return d.parent ? d.children ? \"node\" : \"node node--leaf\" : \"node node--root\"; })
      .style(\"fill\", function(d) { return color(d.color); })      
      .on(\"click\", function(d) { if (focus !== d) zoom(d), d3.event.stopPropagation(); });

  circle = circle.filter(function(d) {return (d.size >= "$filter") && (d.name != \"\")});
//#      .style(\"stroke-dasharray\", \"5,5\", \"important\");
//#      .style(\"stroke\", \"#000\", \"important\")
//#      .style(\"stroke-width\", 1, \"important\");
  
  var text = svg"$i".selectAll(\"text\")
      .data(nodes)
    .enter().append(\"text\")
      .attr(\"class\", \"label\")
      .style(\"fill-opacity\", function(d) { return d.parent === root ? 1 : 0; })
      .style(\"display\", function(d) { return (d.parent === root && d.size >= "$filter") ? \"inline\" : \"none\"; })
      .text( function(d) { return d.name; });

  text = text.filter(function(d) {return (d.size >= "$filter")});
  
  var node = svg"$i".selectAll(\"circle,text\");

  d3.select(\"#"dv_$i"\")
      .style(\"background\", \"rgb(255,255,255)\")
      .on(\"click\", function() { zoom(root); });

  zoomTo([root.x, root.y, root.r * 2 + margin]);

  function zoom(d) {
    var focus0 = focus; focus = d;


    circle.filter(function(d) {return (d.fl == \""$arg"\" && d.parent === focus && focus !== root); })
//#    .style(\"stroke-dasharray\", \"5,5\");
    .style(\"stroke\", \"#000\")
    .style(\"stroke-width\", 4)
    .style(\"stroke-dasharray\", \"4,4\");

    circle.filter(function(d) {return (d.fl == \""$arg"\" && d === focus && focus !== root); })
//#    .style(\"stroke-dasharray\", \"10,10\");
// #    .style(\"stroke\", \"rgb(255,255,255)\")
    .style(\"stroke-width\", 4)
    .style(\"stroke-dasharray\", \"40,40\");

    circle.filter(function(d) {return (d.fl == \""$arg"\" && (d.parent !== focus && d !== focus)); })
    .style(\"stroke\", \"#000\")
    .style(\"stroke-width\", 1)
    .style(\"stroke-dasharray\", \"none\");
    


    var transition = d3.transition()
        .duration(d3.event.altKey ? 7500 : 750)
        .tween(\"zoom\", function(d) {
          var i = d3.interpolateZoom(view, [focus.x, focus.y, focus.r * 2 + margin]);
          return function(t) { zoomTo(i(t)); };
        });
  

    transition.selectAll(\"text\")
      .filter(function(d) { return (d.size >= "$filter" && d.fl == \""$arg"\" &&
      (d.parent === focus || (d === focus && d.children == null) || this.style.display === \"inline\")); })
        .style(\"fill-opacity\", function(d) { return (d.parent === focus || (d === focus && d.children == null))? 1 : 0; })
        .each(\"start\", function(d) { if (d.parent === focus || (d === focus && d.children == null)) this.style.display = \"inline\"; })
        .each(\"end\", function(d) { if (d.parent !== focus && (d !== focus && d.children == null)) this.style.display = \"none\"; });
  }
  function zoomTo(v) {
    var k = diameter / v[2]; view = v;
    node.attr(\"transform\", function(d) { return \"translate(\" + (d.x - v[0]) * k + \",\" + (d.y - v[1]) * k + \")\"; });
    circle.attr(\"r\", function(d) { return d.r * k; });
  }
});" >>$tmp
i=$((i+1));
done

echo "d3.select(self.frameElement).style(\"height\", diameter + \"px\");

</script>" >>$tmp

    
    
    
    
    
    
open -a Safari $tmp

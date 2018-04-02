var D3Component = (function () {
    function D3Component(div) {
        function isHelixInteraction(helix1, helix2, interaction) {
            var start1 = helix1.start;
            var start2 = helix2.start;
            var end1 = helix1.end;
            var end2 = helix2.end;
            var partner1 = interaction.partner1;
            var partner2 = interaction.partner2;

            return (partner1 >= start1 && partner1 <= end1 && partner2 >= start2 && partner2 <= end2) ||
                (partner1 >= start2 && partner1 <= end2 && partner2 >= start1 && partner2 <= end1);
        }

        this.topologyGraph = function(model) {
            var graph = { nodes : [], links : [] };
            model.tm.forEach(function(tm) {
                graph.nodes.push({ name : 'TM' + (graph.nodes.length + 1), start : tm.start, end : tm.end });
            });

            var interactionTypes = ['halogenBonds', 'hydrogenBonds', 'metalComplexes', 'piCationInteractions',
                'piStackings', 'saltBridges', 'waterBridges'];
            for(var i = 0; i < graph.nodes.length - 1; i++) {
                for(var j = i + 1; j < graph.nodes.length; j++) {
                    var helicesInteract = false;
                    var helix1 = graph.nodes[i];
                    var helix2 = graph.nodes[j];

                    interactionTypes.forEach(function(interactionType) {
                        model[interactionType].forEach(function(interaction) {
                            if(isHelixInteraction(helix1, helix2, interaction)) {
                                helicesInteract = true;
                            }
                        });
                    });

                    if(helicesInteract) {
                        graph.links.push({ source : i, target : j, value : 5 });
                    }
                }
            }

            var width = d3.select(div)[0][0].clientWidth, height = 100, radius = 12;
            var color = d3.scale.category20();
            var force = d3.layout.force()
                .gravity(.05)
                .charge(-800)
                .linkDistance(50)
                .size([width, height]);

            var svg = d3.select(div).append("svg")
                .attr("width", width)
                .attr("height", height);

            force.nodes(graph.nodes)
                .links(graph.links)
                .start();

            var link = svg.selectAll(".link")
                .data(graph.links)
                .enter()
                .append("line")
                .attr("class", "link")
                .style("stroke-width", function(d) { return Math.sqrt(d.value) })
                .style("stroke", function(d) { return d3.rgb('#666') });

            var node = svg.selectAll(".node")
                .data(graph.nodes)
                .enter().append("circle")
                .attr("r", radius - .75)
                .style("fill", function(d) { return color('#51b2e9'); })
                .style("stroke", function(d) { return d3.rgb('#51b2e9').darker() })
                .call(force.drag);

            node.append("title")
                .text(function(d) { return d.name; });

            force
                .nodes(graph.nodes)
                .links(graph.links)
                .on("tick", tick)
                .start();

            function tick() {
                node.attr("cx", function(d) { return d.x = Math.max(radius, Math.min(width - radius, d.x)); })
                    .attr("cy", function(d) { return d.y = Math.max(radius, Math.min(height - radius, d.y)); });

                link.attr("x1", function(d) { return d.source.x; })
                    .attr("y1", function(d) { return d.source.y; })
                    .attr("x2", function(d) { return d.target.x; })
                    .attr("y2", function(d) { return d.target.y; });
            }
        }
    }

    return D3Component;
})();
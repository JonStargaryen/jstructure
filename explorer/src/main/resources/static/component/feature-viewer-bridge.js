var FeatureViewerComponent = (function () {
    var fv;

    function FeatureViewerComponent(div, model, featureSelectedCallback, zoomCallback) {
        fv = new FeatureViewer(model.sequence,
            div,
            {
                showAxis: true,
                showSequence: true,
                brushActive: true,
                toolbar: true,
                bubbleHelp: false,
                zoomMax: 5
            });

        if(model.active.length > 0) {
            var activeSites = [];
            model.active.forEach(function (go) {
                activeSites.push({
                    x: +go.position,
                    y: +go.position,
                    description: go.description
                });
            });
            fv.addFeature({
                data: activeSites,
                name: 'Active Sites',
                className: 'tba1',
                type: 'rect',
                color: '#52b1e9'
            });
        }

        if(model.disulfide.length > 0) {
            var bonds = [];
            model.disulfide.forEach(function(go) {
                bonds.push({
                    x : go.start,
                    y : go.end
                });
            });
            fv.addFeature({
                data : bonds,
                name : 'Disulfide Bond',
                className : 'tba2',
                type : 'path',
                color: '#52b1e9'
            });
        }

        if(model.ptm.length > 0) {
            var ptm = [];
            model.ptm.forEach(function (go) {
                ptm.push({
                    x: go.position,
                    y: go.position,
                    description: go.description,
                    // color : go.type
                })
            });
            fv.addFeature({
                data: ptm,
                name: 'Modification',
                className: 'tba3',
                type: 'rect',
                color: '#52b1e9'
            });
        }

        if(model.mutagenesis.length > 0) {
            var mutations = [];
            model.mutagenesis.forEach(function (go) {
                mutations.push({
                    x: +go.position[0],
                    y: +go.position[0],
                    description: go.original + ' \u2192 ' + go.variation  + '\n' + go.description
                });
            });
            fv.addFeature({
                data: mutations,
                name: 'Mutagenesis Site',
                className: 'tba10',
                type: 'rect',
                color: '#52b1e9'
            });
        }

        if(model.variants.length > 0) {
            var variants = [];
            model.variants.forEach(function (go) {
                variants.push({
                    x: go.position,
                    y: go.position,
                    description: go.description
                });
            });
            fv.addFeature({
                data: variants,
                name: 'Natural Variants',
                className: 'tba8',
                type: 'rect',
                color: '#52b1e9'
            });
        }

        // if(model.groups.length > 0) {
        //     var rasa = [];
        //     var present = false;
        //     model.groups.forEach(function (go) {
        //         if(go.rasa !== 0) {
        //             present = true;
        //         }
        //         rasa.push({
        //             x: rasa.length,
        //             y: go.rasa
        //         });
        //     });
        //     if(present) {
        //         fv.addFeature({
        //             data: rasa,
        //             name: 'Accessible Surface',
        //             className: 'tba4',
        //             type: 'line',
        //             color: '#52b1e9',
        //             interpolation: 'basis'
        //         });
        //     }
        // }

        if(model.sse.length > 0) {
            var sse = [];
            model.sse.forEach(function (go) {
                sse.push({
                    x: go.start,
                    y: go.end,
                    description: go.type.toLowerCase()
                });
            });
            fv.addFeature({
                data: sse,
                name: 'Secondary Structure',
                className: 'tba5',
                type: 'rect',
                color: '#52b1e9'
            });
        }

        if(model.tm.length > 0) {
            var tm = [];
            model.tm.forEach(function(go) {
                tm.push({
                    x : go.start,
                    y : go.end,
                    // description : go.type.toLowerCase()
                    description : 'transmembrane'
                })
            });
            fv.addFeature({
                data : tm,
                name : 'Topology',
                className : 'tba6',
                type : 'rect',
                color: '#52b1e9'
            })
        }

        if(featureSelectedCallback !== 'undefined')
            fv.onFeatureSelected(featureSelectedCallback);

        if(zoomCallback !== 'undefined')
            fv.onZoom(zoomCallback);
    }

    return FeatureViewerComponent;
})();
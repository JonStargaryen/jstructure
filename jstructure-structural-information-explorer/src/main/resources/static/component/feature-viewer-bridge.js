var FeatureViewerComponent = (function () {
    var fv;

    function FeatureViewerComponent(div, model, hoverCallback) {
        fv = new FeatureViewer(model.sequence,
            div,
            {
                showAxis: true,
                showSequence: true,
                brushActive: true,
                toolbar: true,
                bubbleHelp: false,
                zoomMax: 5
            },
            hoverCallback);

        var sse = [];
        model.secondaryStructureElements.forEach(function (go) {
            sse.push({
                x: go.start,
                y: go.end,
                description: go.type.toUpperCase()
            });
        });
        fv.addFeature({
            data: sse,
            name: 'Secondary Structure',
            className: 'tba5',
            type: 'rect',
            color: '#52b1e9'
        });

        var early = [];
        model.earlyResidueNumbers.forEach(function (go) {
            early.push({
                x: +go,
                y: +go,
                description: 'Early Folding'
            });
        });
        fv.addFeature({
            data: early,
            name: 'Early Folding',
            className: 'tba1',
            type: 'rect',
            color: '#52b1e9'
        });

        var activeSites = [];
        model.functionalResidueNumbers.forEach(function (go) {
            activeSites.push({
                x: +go,
                y: +go,
                description: 'Functional'
            });
        });
        fv.addFeature({
            data: activeSites,
            name: 'Active Sites',
            className: 'tba1',
            type: 'rect',
            color: '#52b1e9'
        });

        var avgRmsd = [];
        model.residues.forEach(function (go) {
            avgRmsd.push({
                x: avgRmsd.length,
                y: go.averageRmsdIncrease
            });
        });
        fv.addFeature({
            data: avgRmsd,
            name: 'avg. RMSD',
            className: 'tba4',
            type: 'line',
            color: '#52b1e9',
            interpolation: 'basis'
        });

        var maxRmsd = [];
        model.residues.forEach(function (go) {
            maxRmsd.push({
                x: maxRmsd.length,
                y: go.maximumRmsdIncrease
            });
        });
        fv.addFeature({
            data: maxRmsd,
            name: 'max. RMSD',
            className: 'tba4',
            type: 'line',
            color: '#52b1e9',
            interpolation: 'basis'
        });
    }

    return FeatureViewerComponent;
})();
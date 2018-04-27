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
                x: avgRmsd.length + 1,
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

        var sumRmsd = [];
        model.residues.forEach(function (go) {
            sumRmsd.push({
                x: sumRmsd.length + 1,
                y: go.sumRmsdIncrease
            });
        });
        fv.addFeature({
            data: sumRmsd,
            name: 'sum RMSD',
            className: 'tba4',
            type: 'line',
            color: '#52b1e9',
            interpolation: 'basis'
        });

        // var maxRmsd = [];
        // model.residues.forEach(function (go) {
        //     maxRmsd.push({
        //         x: maxRmsd.length + 1,
        //         y: go.maximumRmsdIncrease
        //     });
        // });
        // fv.addFeature({
        //     data: maxRmsd,
        //     name: 'max. RMSD',
        //     className: 'tba4',
        //     type: 'line',
        //     color: '#52b1e9',
        //     interpolation: 'basis'
        // });

        // var eccount = [];
        // model.residues.forEach(function (go) {
        //     eccount.push({
        //         x: eccount.length + 1,
        //         y: go.eccount
        //     });
        // });
        // fv.addFeature({
        //     data: eccount,
        //     name: 'EC count',
        //     className: 'tba4',
        //     type: 'line',
        //     color: '#52b1e9',
        //     interpolation: 'basis'
        // });

        var cumstrength = [];
        model.residues.forEach(function (go) {
            cumstrength.push({
                x: cumstrength.length + 1,
                y: go.cumstrength
            });
        });
        fv.addFeature({
            data: cumstrength,
            name: 'cumulative strength',
            className: 'tba4',
            type: 'line',
            color: '#52b1e9',
            interpolation: 'basis'
        });

        var conservation = [];
        model.residues.forEach(function (go) {
            conservation.push({
                x: conservation.length + 1,
                y: go.conservation
            });
        });
        fv.addFeature({
            data: conservation,
            name: 'conservation',
            className: 'tba4',
            type: 'line',
            color: '#52b1e9',
            interpolation: 'basis'
        });

        // var avgZ = [];
        // model.residues.forEach(function (go) {
        //     avgZ.push({
        //         x: avgZ.length + 1,
        //         y: go.averageRmsdIncreaseZScore
        //     });
        // });
        // fv.addFeature({
        //     data: avgZ,
        //     name: 'average Z',
        //     className: 'tba4',
        //     type: 'line',
        //     color: '#52b1e9',
        //     interpolation: 'basis'
        // });
        //
        // var maxZ = [];
        // model.residues.forEach(function (go) {
        //     maxZ.push({
        //         x: maxZ.length + 1,
        //         y: go.maximumRmsdIncreaseZScore
        //     });
        // });
        // fv.addFeature({
        //     data: maxZ,
        //     name: 'maxiumum Z',
        //     className: 'tba4',
        //     type: 'line',
        //     color: '#52b1e9',
        //     interpolation: 'basis'
        // });
        //
        // var top = [];
        // model.residues.forEach(function (go) {
        //     top.push({
        //         x: top.length + 1,
        //         y: go.fractionOfTopScoringContacts
        //     });
        // });
        // fv.addFeature({
        //     data: top,
        //     name: 'frac top scoring',
        //     className: 'tba4',
        //     type: 'line',
        //     color: '#52b1e9',
        //     interpolation: 'basis'
        // });
    }

    return FeatureViewerComponent;
})();
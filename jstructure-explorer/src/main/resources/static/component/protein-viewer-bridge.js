var ProteinViewerComponent = (function () {
    var sphereAlpha = 1.0;
    var defaultOptions = {
        width : 400,
        height : 400,
        antialias : true,
        fog : false,
        quality : 'high',

        background : '#313a41',
        white : '#FFFFFF',
        lightGray : '#d1e2e6',
        selectionColor : '#ff0000',
        interactions : {
            sphereSize : 1,
            halogenBonds : [0.247, 1.0, 0.749, sphereAlpha],
            hydrogenBonds : [0.0, 0.0, 1.0, sphereAlpha],
            metalComplexes : [0.031, 0.765, 0.976, sphereAlpha],
            piCationInteractions : [0.059, 0.969, 0.941, sphereAlpha],
            piStackings : [0.031, 0.796, 0.149, sphereAlpha],
            saltBridges : [1.0, 1.0, 0.0, sphereAlpha],
            waterBridges : [0.749, 0.749, 1.0, sphereAlpha]
        }
    };

    var viewer, options;

    function ProteinViewerComponent(div) {
        options = defaultOptions;
        viewer = pv.Viewer(document.getElementById(div.substring(1)), options);

        this.loadStructure = function(model) {
            // store all associated information in a custom container - TODO how efficient is this?
            model.pv = {};
            model.pv.structure = io.pdb(model.pdb);
            mol.assignHelixSheet(model.pv.structure);

            // the protein itself
            model.pv.protein = viewer.cartoon('protein_' + model.id,
                model.pv.structure.select('protein'),
                { color: pv.color.uniform(options.white) }
            );

            // potential ligands
            model.pv.ligands = viewer.ballsAndSticks('ligands_' + model.id,
                model.pv.structure.select('ligand'),
                { color: pv.color.uniform(options.white) }
            );

            // interactions
            var interactionTypes = ['halogenBonds', 'hydrogenBonds', 'metalComplexes', 'piCationInteractions',
                'piStackings', 'saltBridges', 'waterBridges'];
            model.pv.interactions = {};
            interactionTypes.forEach(function(entry) {
                renderInteraction(model, entry);
            });

            viewer.centerOn(model.pv.structure);
            viewer.fitTo(model.pv.structure);
        };

        function renderInteraction(model, name) {
            model[name].forEach(function(entry) {
                // skip local interactions
                if(!entry.longRange)
                    return;
                // if the given type of interactions is present, set a flag for easy retrieval
                model.pv.interactions[name] = true;
                viewer.customMesh(name).addSphere(entry.center,
                    options.interactions.sphereSize,
                    { color : options.interactions[name] }
                );
            });
        }

        this.toggleInteractions = function() {
            //TODO impl
        };

        this.toggleLigands = function() {
            //TODO impl
        };
    }

    return ProteinViewerComponent;
})();
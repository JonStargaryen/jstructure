var ProteinViewerComponent = (function () {
    var defaultOptions = {
        width : 400,
        height : 400,
        antialias : true,
        fog : false,
        quality : 'high',
        background : '#313a41',
        white : '#FFFFFF',
        lightGray : '#d1e2e6',
        selectionColor : '#ff0000'
    };

    var viewer, options, structure, rendering, contact;

    function ProteinViewerComponent(div) {
        options = defaultOptions;

        this.loadStructure = function(model) {
            var element = document.getElementById(div.substring(1));
            element.innerHTML = '';
            viewer = pv.Viewer(element, options);
            structure = io.pdb(model.pdbRepresentation);
            mol.assignHelixSheet(structure);

            rendering = viewer.cartoon('protein_' + model.id,
                structure.select('protein'),
                { color: pv.color.uniform(options.white) });

            viewer.centerOn(structure);
            viewer.fitTo(structure);
        };

        this.colorStructure = function(model, feature) {
            model.scale_min = 0;
            model.scale_max = 0;
            model.residues.forEach(function(r) {
                var value = r[feature];
                if(value < model.scale_min) {
                    model.scale_min = value;
                }
                if(value > model.scale_max) {
                    model.scale_max = value;
                }
            });

            // the protein itself
            rendering.colorBy(new pv.color.ColorOp(function(atom, out, index) {
                var residue = atom.residue();
                var value = 0.0;
                model.residues.forEach(function(r) {
                    if(r.residueIdentifier === residue.num()) {
                        value = r[feature];
                    }
                });

                if(value === 0.0) {
                    out[index + 0] = 1.0;
                    out[index + 1] = 1.0;
                    out[index + 2] = 1.0;
                    out[index + 3] = 1.0;
                } else {
                    // index + 0, index + 1 etc. are the positions in the output array
                    // at which the red (+0), green (+1), blue (+2) and  alpha (+3)
                    // components are to be written.
                    if(value < 0) {
                        c = (1 - (model.scale_min - value) / model.scale_min);
                        out[index + 0] = 1.0;
                        out[index + 1] = 1.0 - c;
                        out[index + 2] = 1.0 - c;
                    } else {
                        c = (1 - (model.scale_max - value) / model.scale_max);
                        out[index + 0] = 1.0 - c;
                        out[index + 1] = 1.0;
                        out[index + 2] = 1.0 - c;
                    }
                    out[index + 3] = 1.0;
                }
            }));

            viewer.requestRedraw();
        };

        this.select = function(res1, res2) {
            if(isContact(res1, res2)) {
                contact = viewer.customMesh("contact");
                var c1 = rendering.select({ rnum : res1 }).center();
                var c2 = rendering.select({ rnum : res2 }).center();
                contact.addTube(c1, c2, 0.5, { cap : true, color : 'red' });
            } else {
                rendering.setSelection(rendering.select({ rnum : res1 }));
            }
            viewer.requestRedraw();
        };

        this.selectByIndex = function(i) {
            rendering.setSelection(rendering.select({ rindices : [i] }));
            viewer.requestRedraw();
        };

        this.deselect = function(res1, res2) {
            if(isContact(res1, res2)) {
                viewer.rm("contact");
            } else {
                var ev = structure.full().createEmptyView();
                rendering.setSelection(ev);
            }
            viewer.requestRedraw();
        };

        function isContact(res1, res2) {
            return res1 && res2;
        }
    }

    return ProteinViewerComponent;
})();
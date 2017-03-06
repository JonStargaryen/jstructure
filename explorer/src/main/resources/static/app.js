'use strict';

(function () {
    var MODULE = angular.module('mpe', ['ngRoute', 'ngResource', 'nvd3', 'ngDropdowns']);

    MODULE.constant('design', {
        defaultColor : '#454b52',
        lighterColor : '#d1e2e6',
        highlightColor : '#52b1e9'
    });

    MODULE.constant('features', [
        { name : 'no coloring', tag : 'none', discrete : true },
        { name : 'relative accessible surface area', tag : 'rasa', discrete : false },
        { name : 'secondary structure element', tag: 'sse', discrete : true }
    ]);

    MODULE.constant('renderModes',  [
        'cartoon', 'ballsAndSticks', 'sline', 'lines', 'trace', 'lineTrace', 'tube', 'spheres'
    ]);

    MODULE.constant('secondaryStructures', [
        'coil', 'bend', 'turn', '\u03C0-helix', '3-10-helix', 'bridge', 'extended', '\u03B1-helix'
    ]);

    MODULE.constant('interactions', {
        'halogenBonds' : '#3FFFBF',
        'hydrogenBonds' : '#0000FF',
        'hydrophobicInteractions' : '#7F7F7F',
        'metalComplexes' : '#8C3F99',
        'piCationInteractions' : '#FF7F00',
        'piStackings' : '#8CB266',
        'saltBridges' : '#FFFF00',
        'waterBridges' : '#BFBFFF' });

    /**
     * navigation
     */
    MODULE.config(['$locationProvider', '$routeProvider', function ($locationProvider, $routeProvider) {
        // hide index.html# in browser url
        $locationProvider.html5Mode({
            enabled: true
        });

        // navigation rules
        $routeProvider.when('/home', {
            templateUrl: '/home.htm',
            controller: 'HomeController'
        });

        $routeProvider.when('/about', {
            templateUrl: '/about.htm'
        });

        $routeProvider.when('/protein/:pdbid', {
            templateUrl: '/protein.htm',
            controller: 'ProteinController'
        });

        $routeProvider.otherwise('/home');
    }]);

    /**
     * on app start
     */
    MODULE.run(['$rootScope', '$http', '$location', 'renderModes', 'features', '$filter', 'secondaryStructures',
        function ($rootScope, $http, $location, renderModes, features, $filter, secondaryStructures) {
        $rootScope.alerts = [];

        // message/alert container
        $rootScope.closeAlert = function (index) {
            $rootScope.alerts.splice(index, 1);
        };

        $rootScope.page = function () {
            return $location.path();
        };

        $rootScope.renderModes = [];
        renderModes.forEach(function(entry) {
            $rootScope.renderModes.push({ text: $filter('splitAtUppercase')(entry),
                index: $rootScope.renderModes.length, raw: entry });
        });

        $rootScope.features = features;
        $rootScope.secondaryStructures = secondaryStructures;
    }]);

    /**
     * handles multiple tabs in the UI
     */
    MODULE.controller('TabController', function () {
        this.tab = 1;

        this.selectTab = function (setTab){
            this.tab = setTab;
        };

        this.isSelected = function (checkTab){
            return this.tab === checkTab;
        };
    });

    MODULE.controller('ProteinController', ['$scope', '$rootScope', '$routeParams', '$http', 'design', 'interactions',
        function($scope, $rootScope, $routeParams, $http, design, interactions) {
        $scope.loading = true;
        $scope.options = {
            renderMode : $scope.renderModes[0],
            coloringFeature : $scope.features[0],
            renderLigands : true,
            renderKinks : false,
            renderMembrane : false
        };
        var initializedInteractions = [];
        var viewer, structure, aminoAcids, ligands, kinks;

        $http.get('/api/proteins/pdbid/' + $routeParams.pdbid).then(function(data) {
            $scope.protein = data.data;
            $rootScope.context = {};
            $rootScope.context.name = $routeParams.pdbid;
            $rootScope.context.title = $scope.protein.title;

            /* compose PDB representation of the protein */
            var rep = "";
            $scope.protein.chains.forEach(function(chain) {
                chain.groups.forEach(function(group) {
                    group.atoms.forEach(function(atom) {
                        rep += atom.pdb + "\n";
                    });
                });
            });
            $scope.protein.pdb = rep;

            console.log($scope.protein);
            visualizeProteinStructure();

            $scope.loading = false;
        }, function (d) {
            $rootScope.alerts.push({ type: 'danger', msg: 'loading protein ' + $routeParams.pdbid +
            ' from server failed with [' + d.status + '] ' + d.statusText });
        });

        /* PV related functions */
        function visualizeProteinStructure() {
            var options = {
                // div to be selected is not visible at that time
                width: 400,
                height: 400,
                antialias: true,
                quality : 'high',
                background:'#313a41',
                // animateTime : 500,
                selectionColor : '#f00',
                fog : false
            };

            // ensure container is empty
            document.getElementById('protein-visualizer').innerHTML = '';

            viewer = pv.Viewer(document.getElementById('protein-visualizer'), options);
            structure = io.pdb($scope.protein.pdb);
            mol.assignHelixSheet(structure);

            initProtein();
            initLigands();
            if(!$scope.options.renderLigands)
                toggle('ligands', false, true);
            initKinks();
            if(!$scope.options.renderKinks)
                toggle('kinks', false, true);
            initMembrane();
            if(!$scope.options.renderMembrane)
                toggle('membrane', false, true);

            viewer.centerOn(structure);
            viewer.autoZoom();
        }

        function findGroup(protein, chainId, num) {
            var residue = {};
            protein.chains.forEach(function(chain) {
                if(chain.id != chainId)
                    return;
                chain.groups.forEach(function(group) {
                    if(group.resn != num)
                        return;
                    residue = group;
                });
            });
            return residue;
        }

        /* render plain PDB structure with chosen style */
        function initProtein() {
            aminoAcids = viewer.cartoon('protein', structure.select('protein'), { color: colorOp });
        }

        function initLigands() {
            ligands = viewer.ballsAndSticks('ligands', structure.select('ligand'),
                { color: pv.color.uniform(design.lighterColor) });
        }

        function initMembrane() {
            $scope.protein.membrane.forEach(function(entry) {
                viewer.customMesh('membrane').addSphere([entry[0], entry[1], entry[2]], 0.5,
                    { color : [0.8, 0.8, 0.8] });
            });
        }

        function initKinks() {
            var kinkIds = [];
            $scope.protein.kinks.forEach(function(entry) {
                kinkIds.push({ chain: entry.chain, rnum: +entry.kinkPosition});
            });
            kinks = viewer.ballsAndSticks('kinks', structure.residueSelect(function(res) {
                for(var i = 0; i < kinkIds.length; i++) {
                    if(kinkIds[i].chain == res.chain().name() && kinkIds[i].rnum == res.num())
                        return true;
                }
                return false;
            }), { color: pv.color.uniform(design.highlightColor) });
        }

        var colorOp = new pv.color.ColorOp(function(atom, out, index) {
            var rgb;

            if($scope.options.coloringFeature.tag === "none") {
                // no coloring requested
                rgb = pv.color.hex2rgb(design.lighterColor);
            } else {
                // read property and transform to rgb
                var residue = atom.residue();
                var group = findGroup($scope.protein, residue.chain().name(), residue.num());
                rgb = hsl2rgb(group[$scope.options.coloringFeature.tag] / 3);
            }

            out[index] = rgb[0];
            out[index + 1] = rgb[1];
            out[index + 2] = rgb[2];
            out[index + 3] = 1.0;
        });

        /* render PLIP interactions */
        function createInteractions(type) {
            console.log('creating interactions for ' + type);
            $scope.protein[type].forEach(function(entry) {
                var atom1 = structure.atom(entry.a1);
                var atom2 = structure.atom(entry.a2);
                var g = viewer.customMesh(type);
                g.addTube(atom1.pos(), atom2.pos(), 0.15, { cap : false, color : interactions[type] });
            });
            initializedInteractions.push(type);
        }

        /* watchers for front-end options */
        /*$scope.$watch('options.renderMode', function(newVal, oldVal){
            if(newVal === oldVal)
                return;
            console.log("render mode changed from '" + oldVal.name + "' to '" + newVal.name + "'");
            viewer.rm('protein');
            initProtein();
        });*/

        $scope.$watch('options.coloringFeature', function(newVal, oldVal){
            if(newVal === oldVal)
                return;
            console.log("coloring feature changed from '" + oldVal.name + "' to '" + newVal.name + "'");
            aminoAcids.colorBy(colorOp);
            viewer.requestRedraw();
        });

        function toggle(context, newVal, oldVal) {
            if(newVal === oldVal)
                return;
            console.log("render '" + context + "' changed from '" + oldVal + "' to '" + newVal + "'");

            if(newVal)
                viewer.show(context);
            else
                viewer.hide(context);
            viewer.requestRedraw();
        }

        $scope.$watch('options.renderLigands', function(newVal, oldVal) {
            toggle('ligands', newVal, oldVal);
        });

        $scope.$watch('options.renderKinks', function(newVal, oldVal) {
            toggle('kinks', newVal, oldVal);
        });

        $scope.$watch('options.renderMembrane', function(newVal, oldVal) {
            toggle('membrane', newVal, oldVal);
        });

        $scope.$watch('options.renderHalogenBonds', function(newVal, oldVal) {
            toggleInteractions(newVal, oldVal, 'halogenBonds');
        });

        $scope.$watch('options.renderHydrogenBonds', function(newVal, oldVal) {
            toggleInteractions(newVal, oldVal, 'hydrogenBonds');
        });

        $scope.$watch('options.renderHydrophobicInteractions', function(newVal, oldVal) {
            toggleInteractions(newVal, oldVal, 'hydrophobicInteractions');
        });

        $scope.$watch('options.renderMetalComplexes', function(newVal, oldVal) {
            toggleInteractions(newVal, oldVal, 'metalComplexes');
        });

        $scope.$watch('options.renderPiCationInteractions', function(newVal, oldVal) {
            toggleInteractions(newVal, oldVal, 'piCationInteractions');
        });

        $scope.$watch('options.renderPiStackings', function(newVal, oldVal) {
            toggleInteractions(newVal, oldVal, 'piStackings');
        });

        $scope.$watch('options.renderSaltBridges', function(newVal, oldVal) {
            toggleInteractions(newVal, oldVal, 'saltBridges');
        });

        $scope.$watch('options.renderWaterBridges', function(newVal, oldVal) {
            toggleInteractions(newVal, oldVal, 'waterBridges');
        });

        function toggleInteractions(newVal, oldVal, type) {
            if(newVal === oldVal)
                return;
            console.log('render ' + type + ' changed from ' + oldVal + ' to ' + newVal);
            if(newVal) {
                if(initializedInteractions.indexOf(type) < 0)
                    createInteractions(type);
                viewer.show(type);
            } else {
                viewer.hide(type);
            }
            viewer.requestRedraw();
        }

        $scope.highlight = function(chain, resn, endResn) {
            // console.log(chain + " " + resn + " " + endResn);
            var selection = { cname : chain };
            if(resn) {
                if (endResn)
                    selection.rnumRange = [+resn, +endResn];
                else
                    selection.rnum = +resn;
            }
            aminoAcids.setSelection(aminoAcids.select(selection));
            ligands.setSelection(ligands.select(selection));
            kinks.setSelection(kinks.select(selection));
            viewer.requestRedraw();
        };

        $scope.deselect = function() {
            var ev = structure.full().createEmptyView();
            aminoAcids.setSelection(ev);
            ligands.setSelection(ev);
            kinks.setSelection(ev);
            viewer.requestRedraw();
        };

        // http://stackoverflow.com/questions/2353211/hsl-to-rgb-color-conversion
        function hsl2rgb(h){
            var r, g, b;

            var hue2rgb = function hue2rgb(p, q, t){
                if(t < 0) t += 1;
                if(t > 1) t -= 1;
                if(t < 1/6) return p + (q - p) * 6 * t;
                if(t < 1/2) return q;
                if(t < 2/3) return p + (q - p) * (2/3 - t) * 6;
                return p;
            };

            var q = 0.7 + 0.6 - 0.7 * 0.6;
            var p = 2 * 0.7 - q;
            r = hue2rgb(p, q, h + 1/3);
            g = hue2rgb(p, q, h);
            b = hue2rgb(p, q, h - 1/3);

            return [r, g, b];
        }
    }]);

    /**
     * fetch all available data for the home view
     */
    MODULE.controller('HomeController', ['$scope', '$location', 'ProteinService',
        function($scope, $location, ProteinService) {
        $scope.allProteins = ProteinService.allProteins;
        $scope.nonRedundantAlphaHelicalProteins = ProteinService.nonRedundantAlphaHelicalProteins;
    }]);

    MODULE.factory('ProteinService', ['$rootScope', '$http', function($rootScope, $http) {
        /* fetch all known protein ids */
        var allProteins = [];
        var nonRedundantAlphaHelicalProteins = [];

        $http.get('/api/proteins/all').then(function(d) {
            d.data.forEach(function (p) {
                allProteins.push(p);
            });
        }, function (d) {
            $rootScope.alerts.push({ type: 'danger', msg: 'loading all proteins from server failed with [' + d.status +
            '] ' + d.statusText });
        });

        $http.get('/api/proteins/alpha_nr').then(function(d) {
            d.data.forEach(function (p) {
                nonRedundantAlphaHelicalProteins.push(p);
            });
        }, function (d) {
            $rootScope.alerts.push({ type: 'danger', msg:
            'loading non-redundant alpha helical proteins from server failed with [' + d.status + '] ' + d.statusText });
        });

        return { "allProteins" : allProteins,
            "nonRedundantAlphaHelicalProteins" : nonRedundantAlphaHelicalProteins };
    }]);

    /* format PV render modes to something nice */
    MODULE.filter('splitAtUppercase', function() {
        return function(input) {
            if(input) {
                return input
                    // insert a space before all caps
                    .replace(/([A-Z])/g, ' $1')
                    .toLowerCase();
            }
        }
    });

    MODULE.filter('formatAuthors', function() {
        return function(input) {
            if(input) {
                if(input.length == 0)
                    return '';
                if(input.length < 2)
                    return omitFullStops(input) + '.';
                return omitFullStops(input[0]) + ' et al.';
            }
        };

        function omitFullStops(input) {
            return input.split('.').join('');
        }
    });

    MODULE.filter('formatPositions', function() {
        return function (input) {
            if (input.length == 1)
                return input[0];
            if (input.length == 2)
                return input[0] + '-' + input[1];
        };
    });
})();
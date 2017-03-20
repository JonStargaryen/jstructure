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
        // { name : 'amino acid', tag : 'aa', discrete : true },
        // { name : 'relative accessible surface area', tag : 'rasa', discrete : false },
        // { name : 'secondary structure element', tag : 'sse', discrete : true },
        // { name : 'spatial proximity', tag : 'cerosene', discrete : false }

    ]);

    MODULE.constant('renderModes',  [
        'cartoon', 'ballsAndSticks', 'sline', 'lines', 'trace', 'lineTrace', 'tube', 'spheres'
    ]);

    MODULE.constant('secondaryStructures', [
        'coil', 'bend', 'turn', '\u03C0-helix', '3-10-helix', 'bridge', 'extended', '\u03B1-helix'
    ]);

    MODULE.constant('aminoAcids', [
        { name : 'alanine', olc : 'A', tlc : 'Ala' },
        { name : 'arginine', olc : 'R', tlc : 'Arg' },
        { name : 'asparagine', olc : 'N', tlc : 'Asn' },
        { name : 'aspartic acid', olc : 'D', tlc : 'Asp' },
        { name : 'cysteine', olc : 'C', tlc : 'Cys' },
        { name : 'glutamic acid', olc : 'E', tlc : 'Glu' },
        { name : 'glutamine', olc : 'Q', tlc : 'Gln' },
        { name : 'glycine', olc : 'G', tlc : 'Gly' },
        { name : 'histidine', olc : 'H', tlc : 'His' },
        { name : 'isoleucine', olc : 'I', tlc : 'Ile' },
        { name : 'leucine', olc : 'L', tlc : 'Leu' },
        { name : 'lysine', olc : 'K', tlc : 'Lys' },
        { name : 'methionine', olc : 'M', tlc : 'Met' },
        { name : 'phenylalanine', olc : 'F', tlc : 'Phe' },
        { name : 'proline', olc : 'P', tlc : 'Pro' },
        { name : 'serine', olc : 'S', tlc : 'Ser' },
        { name : 'threonine', olc : 'T', tlc : 'Thr' },
        { name : 'tryptophan', olc : 'W', tlc : 'Trp' },
        { name : 'tyrosine', olc : 'Y', tlc : 'Tyr' },
        { name : 'valine', olc : 'V', tlc : 'Val' }
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

        $routeProvider.when('/chain/:chainid', {
            templateUrl: '/protein.htm',
            controller: 'ProteinController'
        });

        $routeProvider.when('/literature', {
            templateUrl: '/literature.htm'
        });

        $routeProvider.otherwise('/home');
    }]);

    /**
     * on app start
     */
    MODULE.run(['$rootScope', '$http', '$location', 'renderModes', 'features', '$filter', 'secondaryStructures', 'aminoAcids',
        function ($rootScope, $http, $location, renderModes, features, $filter, secondaryStructures, aminoAcids) {
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
        $rootScope.aminoAcids = aminoAcids;
    }]);

    MODULE.controller('MultiSequenceController', ['$scope', 'ProteinService',
        function($scope, ProteinService) {
        var initialized = false;
        $scope.sequences = [];

        $scope.$watch('tab', function(newVal) {
            // tab got id 1
            if(!initialized &&newVal === 1) {
                initialize();
                initialized = true;
                console.log($scope.sequences);
            }
        });

        function initialize() {
            /* load associated sequences from server */
            ProteinService.loadSequences($scope.protein.rep).then(function(response) {
                $scope.sequences = response.data.sequences;

                var MyOptions = Array();
                MyOptions.border = false;
                MyOptions.nocolor = true;
                MyOptions.scrollX = "100%";
                MyOptions.colourScheme = "zappo";
                MyOptions.selectable = false;
                MyOptions.plainTooltips = true;
                printJSAV('sequenceDisplay', $scope.sequences, MyOptions);
                console.log($scope.sequences);
            }).catch(function(response) {
                console.log('impl error handling at multi-sequence view:');
                console.log(response);
            });
        }
    }]);

    MODULE.controller('ProteinController', ['$scope', '$rootScope', '$routeParams', 'design', 'interactions', 'ProteinService',
        function($scope, $rootScope, $routeParams, design, interactions, ProteinService) {
        $scope.loading = true;
        $scope.options = {
            renderMode : $scope.renderModes[0],
            coloringFeature : $scope.features[0],
            renderLigands : true
        };
        var initializedInteractions = [];
        var viewer, structure, protein, ligands;

        $scope.tab = 0;
        $scope.selectTab = function (setTab){
            $scope.tab = setTab;
        };
        $scope.isSelected = function (checkTab){
            return $scope.tab === checkTab;
        };

        ProteinService.loadModel($routeParams.chainid).then(function(response) {
            $scope.protein = response.data;
            console.log($scope.protein);

            // the context navigation in the header - push entries for additional levels
            $rootScope.context = [];
            var split = $routeParams.chainid.split("_");
            $scope.protein.pdbId = split[0];
            $scope.protein.chainId = split[1];
            $rootScope.context.push('pdb: ' + $scope.protein.pdbId);
            $rootScope.context.push('chain: ' + $scope.protein.chainId);

            visualizeProteinStructure();

            // $scope.protein.sequences = [];
            // $scope.protein.homologous.forEach(function(homolog) {
            //     ProteinService.loadSequence(homolog).then(function(response) {
            //         $scope.protein.sequences.push(response.data);
            //     }).catch(function(response) {
            //         alert('sequences down with ' + response);
            //     });
            // });
            //
            // console.log($scope.protein.sequences);

            // var MyOptions = Array();
            // MyOptions.border = false;
            // MyOptions.highlight = [3,5,10,14];
            // MyOptions.fasta = true;
            // MyOptions.consensus = true;
            // MyOptions.colourScheme = "zappo";
            // MyOptions.plainTooltips = true;
            // printJSAV('sequenceDisplay', $scope.protein.sequences, MyOptions);

            $scope.loading = false;
        }).catch(function(response) {
            alert('down');
        });

        /* PV related functions */
        function visualizeProteinStructure() {
            var options = {
                width : 400,
                height : 400,
                // antialias : true,
                quality : 'high',
                background : '#313a41',
                animateTime : 500,
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
            protein = viewer.cartoon('protein', structure.select('protein'), { color: colorOp });
        }

        function initLigands() {
            ligands = viewer.ballsAndSticks('ligands', structure.select('ligand'),
                { color: pv.color.uniform(design.lighterColor) });
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
                var hsv = group[$scope.options.coloringFeature.tag];
                rgb = hsl2rgb(hsv[0] / 360, hsv[1] / 100, hsv[2] / 100);
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
                var g = viewer.customMesh(type);
                g.addTube(entry.coords1, entry.coords2, 0.15, { cap : false, color : interactions[type] });
            });
            initializedInteractions.push(type);
        }

        $scope.$watch('options.coloringFeature', function(newVal, oldVal){
            if(newVal === oldVal)
                return;
            console.log("coloring feature changed from '" + oldVal.name + "' to '" + newVal.name + "'");
            protein.colorBy(colorOp);
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
            protein.setSelection(protein.select(selection));
            ligands.setSelection(ligands.select(selection));
            viewer.requestRedraw();
        };

        $scope.deselect = function() {
            var ev = structure.full().createEmptyView();
            protein.setSelection(ev);
            ligands.setSelection(ev);
            viewer.requestRedraw();
        };

        // http://stackoverflow.com/questions/2353211/hsl-to-rgb-color-conversion
        function hsl2rgb(h, s, v){
            var r, g, b;

            if(s == 0) {
                r = g = b = v; // achromatic
            } else {
                var hue2rgb = function hue2rgb(p, q, t) {
                    if (t < 0) t += 1;
                    if (t > 1) t -= 1;
                    if (t < 1 / 6) return p + (q - p) * 6 * t;
                    if (t < 1 / 2) return q;
                    if (t < 2 / 3) return p + (q - p) * (2 / 3 - t) * 6;
                    return p;
                };

                var q = v < 0.5 ? v * (1 + s) : v + s - v * s;
                var p = 2 * v - q;
                r = hue2rgb(p, q, h + 1 / 3);
                g = hue2rgb(p, q, h);
                b = hue2rgb(p, q, h - 1 / 3);
            }

            return [r, g, b];
        }
    }]);

    /**
     * fetch all available data for the home view
     */
    MODULE.controller('HomeController', ['$scope', '$location', 'ProteinService',
        function($scope, $location, ProteinService) {
        $scope.allChains = ProteinService.allChains;
        $scope.representativeChains = ProteinService.representativeChains;
    }]);

    MODULE.factory('ProteinService', ['$rootScope', '$http', function($rootScope, $http) {
        var apiUrl = '/api/chains/';

        var reps = [];
        var all = [];

        function logError(context, error) {
            $rootScope.alerts.push({ type: 'danger',
                msg: 'loading ' + context + ' from server failed with [' + error.status + '] ' + error.statusText
            });
        }

        $http.get(apiUrl + 'reps').then(function(response) {
            response.data.forEach(function(id) {
                reps.push(id);
            });
        }, function(response) {
            logError('representative chain ids', response)
        });
        $http.get(apiUrl + 'all').then(function(response) {
            response.data.forEach(function(id) {
                all.push(id);
            });
        }, function(response) {
            logError('all chain ids', response)
        });

        return {
            allChains : all,
            representativeChains : reps,
            loadSequences : function(chainId) {
                return $http({ url : apiUrl + 'alignment/' + chainId, method : 'GET', responseType: 'text' });
            },
            loadStructure : function(chainId) {
                return $http({ url : apiUrl + 'structure/' + chainId, method : 'GET', responseType: 'text' });
            },
            loadModel : function(chainId) {
                return $http.get(apiUrl + 'json/' + chainId);
            }
        }
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
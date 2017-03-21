'use strict';

(function () {
    var MODULE = angular.module('mpe', ['ngRoute', 'ngResource', 'nvd3', 'ngDropdowns']);

    MODULE.constant('design', {
        defaultColor : '#454b52',
        lighterColor : '#d1e2e6',
        highlightColor : '#52b1e9'
    });

    MODULE.constant('features', [
        { name : 'no coloring', tag : 'none', discrete : true }
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
    MODULE.run(['$rootScope', '$http', '$location',
        function ($rootScope, $http, $location) {
        $rootScope.alerts = [];

        // message/alert container
        $rootScope.closeAlert = function (index) {
            $rootScope.alerts.splice(index, 1);
        };

        $rootScope.page = function () {
            return $location.path();
        };
    }]);

    MODULE.controller('MultiSequenceController', ['$scope', 'ProteinService',
        function($scope, ProteinService) {
        var dataInitialized = false;
        // positions with variants
        var variants = [];
        // collection of already initialized (i.e. PDB fetched from server) chains
        var initializedChains = [];
        // collection of currently selected chains
        $scope.selectedChains = {};

        $scope.$watch('tab', function(newVal) {
            // tab got id 1 TODO this is ugly
            if(!dataInitialized &&newVal === 1) {
                initialize();
                dataInitialized = true;
            }
        });

        function initialize() {
            // load associated sequences from server
            ProteinService.loadSequences($scope.protein.rep).then(function(response) {
                var alignment = response.data;
                console.log(alignment);

                // expose sequences
                $scope.chains = alignment.chains;

                var sequenceLength = $scope.chains[0].sequence.length;
                var numberOfSequences = $scope.chains.length;
                // determine positions with variants
                outer : for(var pos = 0; pos < sequenceLength; pos++) {
                    var char = $scope.chains[0].sequence[pos].olc;
                    for(var seq = 1; seq < numberOfSequences; seq++) {
                        if(char != $scope.chains[seq].sequence[pos].olc) {
                            //TODO some entropy value here?
                            variants.push(pos);
                            continue outer;
                        }
                    }
                }

                // mark representative as selected
                initializedChains.push(alignment.representativeChainId);
                $scope.selectedChains[alignment.representativeChainId] = true;
            }).catch(function(response) {
                console.log('impl error handling at multi-sequence view:');
                console.log(response);
            });
        }

        $scope.isVariant = function(index) {
            return variants.indexOf(index) != -1;
        };

        $scope.toggleChain = function(index) {
            var chainId = $scope.chains[index].id;
            // ensure chain is loaded
            if(initializedChains.indexOf(chainId) == -1) {
                initializeChain(chainId);
            }

            //TODO impl add/remove, show/hide in pv
        };

        function initializeChain(chainId) {
            console.log("loading " + chainId);

            ProteinService.loadStructure(chainId).then(function(response) {
                console.log(response.data);
            }).catch(function(response) {
                console.log("TODO error handling:");
                console.log(response);
            });

            initializedChains.push(chainId);
        }
    }]);

    MODULE.controller('ProteinController', ['$scope', '$rootScope', '$routeParams', 'design', 'interactions', 'ProteinService', 'features',
        function($scope, $rootScope, $routeParams, design, interactions, ProteinService, features) {
        $scope.loading = true;
        $scope.options = {
            coloringFeature : features[0],
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
            console.log("highlighting: " + chain + " " + resn + " " + endResn);
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

    /**
     * access to the back-end
     */
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
                return $http({ url : apiUrl + 'pdb/' + chainId, method : 'GET', responseType: 'text' });
            },
            loadModel : function(chainId) {
                return $http.get(apiUrl + 'json/' + chainId);
            }
        }
    }]);
})();
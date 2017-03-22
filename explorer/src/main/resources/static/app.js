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
            templateUrl: '/chain.htm',
            controller: 'ChainController'
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

        $rootScope.$on("$locationChangeStart", function() {
           $rootScope.context = [];
        });
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
     * the controller for the chain-view
     */
    MODULE.controller('ChainController', ['$scope', '$rootScope', '$routeParams', 'ProteinService', 'design',
        function($scope, $rootScope, $routeParams, ProteinService, design) {
        // loading flags
        $scope.loadingModel = true;
        $scope.loadingAlignment = true;
        $scope.loadingStructure = true;

        // the full-fledged model instance
        $scope.reference = {};

        // map of pdb structures: id -> pdb-record
        $scope.chains = {};

        // the sequence alignment and variant positions therein
        $scope.alignment = {};
        $scope.variants = [];
        //TODO impl
        $scope.annotations = [0, 20, 30];
        var sequenceLength = 0;
        var focusPosition = -1000;
        var hoverRange = { start : -1000, end : -1000 };
        $scope.selectionMode = false;
        $scope.windowSize = 3;

        // instance and control over pv
        var viewerDefaultOptions = {
            width : 400,
            height : 400,
            antialias : true,
            quality : 'high',
            background : '#313a41',
            selectionColor : '#f00',
            fog : false
        };
        var viewer = pv.Viewer(document.getElementById('protein-visualizer'), viewerDefaultOptions);
        $scope.viewerOptions = {
            renderLigands : false
        };

        // control over the protein-view navigation tab
        $scope.tab = 0;
        $scope.selectTab = function (setTab) {
            $scope.tab = setTab;
        };
        $scope.isSelected = function (checkTab) {
            return $scope.tab === checkTab;
        };

        // init model for reference chain
        console.log('fetching model for ' + $routeParams.chainid);
        ProteinService.loadModel($routeParams.chainid).then(function (response) {
            var chainId = $routeParams.chainid;
            $scope.reference = response.data;
            console.log($scope.reference);

            // the context navigation in the header - push entries for additional levels
            var split = chainId.split("_");
            $scope.reference.pdbId = split[0];
            $scope.reference.chainId = split[1];

            // register in model map
            $scope.chains[chainId] = $scope.reference;
            $scope.chains[chainId].selected = true;

            registerInViewer(chainId, true);
            initializeAlignment();
        }).catch(function (response) {
            console.log('model: impl error handling!');
            console.log(response);
        }).finally(function () {
            $scope.loadingModel = false;
        });

        function registerInViewer(chainId, center) {
            var chain = $scope.chains[chainId];
            chain.structure = io.pdb(chain.pdb);
            mol.assignHelixSheet(chain.structure);

            chain.pv_protein = viewer.cartoon('protein_' + chainId, chain.structure.select('protein'), { color: pv.color.uniform(design.lighterColor) });
            chain.pv_ligands = viewer.ballsAndSticks('ligands_' + chainId, chain.structure.select('ligand'), { color: pv.color.uniform(design.lighterColor) });
            if(!$scope.viewerOptions.renderLigands) {
                viewer.hide('ligands_' + chainId);
            }

            if(center) {
                viewer.centerOn(chain.structure);
                viewer.autoZoom();
            }
        }

        // fetch additional model instance from the server
        function initalizeStructure(chainId) {
            $scope.loadingStructure = true;
            // load pdb structure
            console.log('fetching additional model for ' + chainId);
            ProteinService.loadModel(chainId).then(function(response) {
                $scope.chains[chainId] = response.data;
                var chain = $scope.chains[chainId];
                chain.selected = true;
                registerInViewer(chainId);
                // duplicated as newly fetched chains are always highlighted
                chain.pv_protein.setSelection(chain.pv_protein.select({ cname : chain.chainId }));
                chain.pv_ligands.setSelection(chain.pv_ligands.select({ cname : chain.chainId }));
            }).catch(function(response) {
                console.log('pdb: impl error handling!');
                console.log(response);
            }).finally(function() {
                $scope.loadingStructure = false;
            });
        }

        function initializeAlignment() {
            // load associated sequences from server
            ProteinService.loadAlignment($scope.reference.rep).then(function(response) {
                $scope.alignment = response.data;

                sequenceLength = $scope.alignment.chains[0].sequence.length;
                var numberOfSequences = $scope.alignment.chains.length;
                // determine positions with variants
                outer : for(var pos = 0; pos < sequenceLength; pos++) {
                    var char = $scope.alignment.chains[0].sequence[pos].olc;
                    for(var seq = 1; seq < numberOfSequences; seq++) {
                        if(char != $scope.alignment.chains[seq].sequence[pos].olc) {
                            $scope.variants.push(pos);
                            continue outer;
                        }
                    }
                }
            }).catch(function(response) {
                console.log('impl error handling at multi-sequence view:');
                console.log(response);
            }).finally(function() {
                $scope.loadingAlignment = false;
            });
        }

        /* pv controls */
        $scope.$watch('viewerOptions.renderLigands', function(newVal, oldVal) {
            if(newVal === oldVal)
                return;

            forAllChains(function(chain) {
                var id = 'ligands_' + chain.id;
                if(newVal && chain.selected) {
                    viewer.show(id);
                } else {
                    viewer.hide(id);
                }
            });
            viewer.requestRedraw();
        });

        /* interactivity functions from the msa-view */

        $scope.selectPosition = function(index) {
            $scope.selectionMode = true;
            // we update potentially once again, to ensure correct ranges in all cases
            focusPosition = index;
            determineHoverRange(index);
            $rootScope.context = ['selection mode', 'range: [' + hoverRange.start + ', ' + hoverRange.end + ']'];
            viewer.hide('*');

            var selection = null;
            forAllChains(function(chain) {
                if(!chain.selected) {
                    return;
                }
                selection = viewer.ballsAndSticks('selection',
                    chain.pv_protein.select({ rnumRange : [ +chain.aminoAcids[hoverRange.start].resn,
                            +chain.aminoAcids[hoverRange.end].resn ] }),
                    { color: pv.color.uniform(design.lighterColor) });
            });

            //TODO bugged! is null
            viewer.centerOn(selection);
            viewer.autoZoom();
            viewer.requestRedraw();
        };

        $scope.leaveSelectionMode = function() {
            $scope.selectionMode = false;
            $rootScope.context = [];
            viewer.rm('selection');
            forAllChains(function(chain) {
                if(!chain.selected) {
                    return;
                }
                viewer.show('protein_' + chain.id);
                if($scope.viewerOptions.renderLigands) {
                    viewer.show('ligands_' + chain.id);
                }
            });
            viewer.requestRedraw();
        };

        /* hover selection of 1 residue, to be rendered as balls-and-sticks */
        $scope.hoverPosition = function(index) {
            // do not visualize anything in selection mode
            if($scope.selectionMode) {
                return;
            }

            focusPosition = index;
            determineHoverRange(index);

            forAllChains(function(chain) {
                if(chain.selected == false) {
                    return;
                }

                //TODO there is a offset for gaps
                //TODO performance: faster when rendering selecting as selected style, then creating new objects?
                // console.log(chain.id + ' index: ' + index + ' is rnum: ' + chain.aminoAcids[index].resn);
                viewer.ballsAndSticks('hover',
                    chain.pv_protein.select({ rnumRange : [ +chain.aminoAcids[hoverRange.start].resn,
                        +chain.aminoAcids[hoverRange.end].resn ] }),
                    { color: pv.color.uniform(design.lighterColor) }
                );
            });
        };

        function determineHoverRange(index) {
            hoverRange.start = index - $scope.windowSize;
            if(hoverRange.start < 0) {
                hoverRange.start = 0;
            }
            hoverRange.end = index + $scope.windowSize;
            if(hoverRange.end > sequenceLength) {
                hoverRange.end = sequenceLength;
            }
        }

        $scope.isFocusPosition = function(index) {
            return focusPosition == index;
        };

        $scope.isHoverPosition = function(index) {
            return index >= hoverRange.start && index <= hoverRange.end;
        };

        $scope.dehoverPosition = function() {
            // ignore when in selection mode
            if($scope.selectionMode) {
                return;
            }
            focusPosition = -1;
            hoverRange.start = -1000;
            hoverRange.end = -1000;
            viewer.rm('hover');
            viewer.requestRedraw();
        };

        // toggle chain in pv, potentially fetch it first
        $scope.toggleChain = function(index) {
            var chainId = $scope.alignment.chains[index].id;
            var chain = $scope.chains[chainId];
            // ensure chain is loaded
            if(!chain || !chain.hasOwnProperty('pdb')) {
                initalizeStructure(chainId);
                return;
            }
            chain.selected = !chain.selected;

            if(chain.selected) {
                viewer.show('protein_' + chainId);
                if($scope.viewerOptions.renderLigands) {
                    viewer.show('ligands_' + chainId);
                }
            } else {
                viewer.hide('protein_' + chainId);
                viewer.hide('ligands_' + chainId);
            }
            viewer.requestRedraw();
        };

        $scope.hoverChain = function(index) {
            var chainId = $scope.alignment.chains[index].id;
            var chain = $scope.chains[chainId];
            if(!chain || !chain.selected) {
                return;
            }
            chain.pv_protein.setSelection(chain.pv_protein.select({ cname : chain.chainId }));
            chain.pv_ligands.setSelection(chain.pv_ligands.select({ cname : chain.chainId }));
            viewer.requestRedraw();
        };

        $scope.dehoverChain = function() {
            forAllChains(function(chain) {
                if(!chain.selected) {
                    return;
                }
                var ev = chain.structure.full().createEmptyView();
                chain.pv_protein.setSelection(ev);
                chain.pv_ligands.setSelection(ev);
            });
            viewer.requestRedraw();
        };

        // convenience for-all-chains function
        function forAllChains(func) {
            for(var chainEntry in $scope.chains) {
                if(!$scope.chains.hasOwnProperty(chainEntry)) {
                    return;
                }
                var chain = $scope.chains[chainEntry];
                func(chain);
            }
        }
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
            loadAlignment : function(chainId) {
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
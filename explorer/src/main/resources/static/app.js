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
    MODULE.controller('ChainController', ['$scope', '$rootScope', '$routeParams', 'ProteinService', 'design', 'features',
        function($scope, $rootScope, $routeParams, ProteinService, design, features) {
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
        $scope.features = features;
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
            chain.pv_structure = io.pdb(chain.pdb);
            mol.assignHelixSheet(chain.pv_structure);

            chain.pv_protein = viewer.cartoon('protein_' + chainId, chain.pv_structure.select('protein'), { color: pv.color.uniform(design.highlightColor) });
            chain.pv_ligands = viewer.ballsAndSticks('ligands_' + chainId, chain.pv_structure.select('ligand'), { color: pv.color.uniform(design.lighterColor) });
            if(!$scope.viewerOptions.renderLigands || $scope.selectionMode) {
                viewer.hide('ligands_' + chainId);
            }
            if($scope.selectionMode) {
                viewer.hide('protein_' + chainId);
                viewer.show('selection_' + chainId);
                // new structure registered while in selection mode: assign pv_selection property
                chain.pv_selection = viewer.ballsAndSticks('selection_' + chain.id,
                    chain.pv_protein.select({ rnumRange : [hoverRange.start, hoverRange.end] }),
                    { color: pv.color.uniform(design.highlightColor) });
            }

            if(center) {
                viewer.centerOn(chain.pv_structure);
                viewer.fitTo(chain.pv_structure);
            }
        }

        // fetch additional model instance from the server
        function initializeStructure(chainId) {
            $scope.loadingStructure = true;
            // load pdb pv_structure
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
                sequenceLength = $scope.alignment.sequences[0].sequence.length;
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

        $scope.enterSelectionMode = function(index) {
            $scope.selectionMode = true;
            // we update potentially once again, to ensure correct ranges in all cases
            focusPosition = index;
            determineHoverRange(index);
            $rootScope.context = ['selection mode', 'range: [' + hoverRange.start + ', ' + hoverRange.end + ']'];
            viewer.hide('*');

            forAllChains(function(chain) {
                viewer.rm('selection_' + chain.id);
                chain.pv_selection = viewer.ballsAndSticks('selection_' + chain.id,
                    chain.pv_protein.select({ rnumRange : [hoverRange.start, hoverRange.end] }),
                    { color: pv.color.uniform(design.highlightColor) });
                if(!chain.selected) {
                    viewer.hide('selection_' + chain.id);
                }
                if($scope.viewerOptions.renderLigands) {
                    viewer.show('ligands_' + chain.id);
                }
            });

            var center = $scope.reference.pv_selection.select({ cname : $scope.reference.chainId });
            viewer.centerOn(center);
            viewer.fitTo(center);
            // some safety net as sometimes the highlight effect lingers after this call returns
            $scope.dehoverChain();
            viewer.requestRedraw();
        };

        $scope.leaveSelectionMode = function() {
            $scope.selectionMode = false;
            $rootScope.context = [];
            forAllChains(function(chain) {
                viewer.rm('selection_' + chain.id);
                if(!chain.selected) {
                    return;
                }
                viewer.show('protein_' + chain.id);
                if($scope.viewerOptions.renderLigands) {
                    viewer.show('ligands_' + chain.id);
                }
            });
            viewer.centerOn($scope.reference.pv_structure.select({ cname: $scope.reference.chainId }));
            viewer.autoZoom();
            // some safety net as sometimes the highlight effect lingers after this call returns
            $scope.dehoverChain();
            viewer.requestRedraw();
        };

        /* hover selection of residues, to be rendered as balls-and-sticks */
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

                //TODO performance: faster when rendering selecting as selected style than creating new objects?
                viewer.ballsAndSticks('hover',
                    chain.pv_protein.select({ rnumRange : [hoverRange.start, hoverRange.end] }),
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
            var chainId = $scope.alignment.sequences[index].id;
            var chain = $scope.chains[chainId];
            // ensure chain is loaded
            if(!chain || !chain.hasOwnProperty('pdb')) {
                initializeStructure(chainId);
                return;
            }
            chain.selected = !chain.selected;

            if(chain.selected) {
                // either show full protein or selected fragment
                viewer.show(($scope.selectionMode ? 'selection_' : 'protein_') + chainId);
                if($scope.viewerOptions.renderLigands) {
                    viewer.show('ligands_' + chainId);
                }
            } else {
                viewer.hide('protein_' + chainId);
                viewer.hide('ligands_' + chainId);
                viewer.hide('selection_' + chainId);
            }
            viewer.requestRedraw();
        };

        $scope.hoverChain = function(index) {
            var chainId = $scope.alignment.sequences[index].id;
            var chain = $scope.chains[chainId];
            if(!chain || !chain.selected) {
                return;
            }
            if($scope.selectionMode) {
                chain.pv_selection.setSelection(chain.pv_selection.select({ cname : chain.chainId }));
            } else {
                chain.pv_protein.setSelection(chain.pv_protein.select({ cname : chain.chainId }));
                chain.pv_ligands.setSelection(chain.pv_ligands.select({ cname : chain.chainId }));
            }
            viewer.requestRedraw();
        };

        $scope.dehoverChain = function() {
            forAllChains(function(chain) {
                if(!chain.selected) {
                    return;
                }
                var ev = chain.pv_structure.full().createEmptyView();
                chain.pv_protein.setSelection(ev);
                chain.pv_ligands.setSelection(ev);
                if(chain.pv_selection != undefined) {
                    chain.pv_selection.setSelection(ev);
                }
            });
            viewer.requestRedraw();
        };

        // convenience for-all-chains function
        function forAllChains(func) {
            for(var chainEntry in $scope.chains) {
                if(!$scope.chains.hasOwnProperty(chainEntry)) {
                    return;
                }
                func($scope.chains[chainEntry]);
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
            loadModel : function(chainId) {
                return $http.get(apiUrl + 'json/' + chainId);
            }
        }
    }]);
})();
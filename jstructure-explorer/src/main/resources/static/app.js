'use strict';

(function () {
    var MODULE = angular.module('mpe', ['ngRoute', 'ngResource', 'ngDropdowns']);

    MODULE.constant('design', {
        defaultColor : '#454b52',
        lighterColor : '#d1e2e6',
        highlightColor : '#52b1e9'
    });

    MODULE.constant('features', [
        { name : 'employ no coloring', tag : 'none', discrete : true }
    ]);

    MODULE.constant('interactionTypes', [
        { name : 'visualize no interactions', tag : 'none' },
        { name : 'halogen bonds', tag : 'halogenBonds' },
        { name : 'hydrogen bonds', tag : 'hydrogenBonds' },
        { name : 'metal complexes', tag : 'metalComplexes' },
        { name : 'pi-cation interactions', tag : 'piCationInteractions' },
        { name : 'pi-stacking interactions', tag : 'piStackings' },
        { name : 'salt bridges', tag : 'saltBridges' },
        { name : 'water bridges', tag : 'waterBridges' }
    ]);

    MODULE.constant('aminoAcids', [
        { olc : 'A', tlc : 'Ala', name : 'alanine', gutteridge : 'none' },
        { olc : 'R', tlc : 'Arg', name : 'arginine', gutteridge : 'guanidinium' },
        { olc : 'N', tlc : 'Asn', name : 'asparagine', gutteridge : 'amide' },
        { olc : 'D', tlc : 'Asp', name : 'aspartic acid', gutteridge : 'carboxylate' },
        { olc : 'C', tlc : 'Cys', name : 'cysteine', gutteridge : 'thiol' },
        { olc : 'Q', tlc : 'Gln', name : 'glutamine', gutteridge : 'amide' },
        { olc : 'E', tlc : 'Glu', name : 'glutamic acid', gutteridge : 'carboxylate' },
        { olc : 'G', tlc : 'Gly', name : 'glycine', gutteridge : 'none' },
        { olc : 'H', tlc : 'His', name : 'histidine', gutteridge : 'imidazole'},
        { olc : 'I', tlc : 'Ile', name : 'isoleucine', gutteridge : 'none' },
        { olc : 'L', tlc : 'Leu', name : 'leucine', gutteridge : 'none' },
        { olc : 'K', tlc : 'Lys', name : 'lysine', gutteridge : 'amino' },
        { olc : 'M', tlc : 'Met', name : 'methionine', gutteridge : 'thiol' },
        { olc : 'F', tlc : 'Phe', name : 'phenylalanine', gutteridge : 'none' },
        { olc : 'P', tlc : 'Pro', name : 'proline', gutteridge : 'none' },
        { olc : 'S', tlc : 'Ser', name : 'serine', gutteridge : 'hydroxyl' },
        { olc : 'T', tlc : 'Thr', name : 'threonine', gutteridge : 'hydroxyl' },
        { olc : 'W', tlc : 'Trp', name : 'tryptophan', gutteridge : 'none' },
        { olc : 'Y', tlc : 'Tyr', name : 'tyrosine', gutteridge : 'hydroxyl' },
        { olc : 'V', tlc : 'Val', name : 'valine', gutteridge : 'none' }
    ]);

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
    MODULE.controller('ChainController', ['$scope', '$rootScope', '$routeParams', 'ProteinService', 'design', 'features', 'interactionTypes', 'aminoAcids',
        function($scope, $rootScope, $routeParams, ProteinService, design, features, interactionTypes, aminoAcids) {
        var proteinViewer, featureViewer;

        // var step = 0;
        // $scope.mutationResidue = {};
        // $scope.alternativeResidues = [];
        // $scope.alternativeResidue = {};
        // $scope.mutationDescription = '';
        //
        // $scope.onStep = function(stepToEvaluate) {
        //     return step === stepToEvaluate;
        // };
        //
        // $scope.moveToMutationStep = function() {
        //     if(!$scope.mutationResidue)
        //         return;
        //
        //     step++;
        //     aminoAcids.forEach(function(go) {
        //         if(go.olc === $scope.mutationResidue.aa)
        //             return;
        //         $scope.alternativeResidues.push(go);
        //     });
        // };
        //
        // $scope.mutate = function() {
        //     step++;
        //     $scope.mutationDescription = "mutating '" + $scope.mutationResidue.aa + "-" + $scope.mutationResidue.resn +
        //         "' to '" + $scope.alternativeResidue.name + "'";
        //     $scope.mutationStatus = 'in progress...';
        //
        //     ProteinService.mutate($scope.reference.id, $scope.mutationResidue.resn, $scope.alternativeResidue.olc).then(function(response) {
        //         var data = response.data;
        //         console.log(data);
        //         $scope.mutation = {
        //             original : data.original,
        //             mutation : data.mutation
        //         };
        //
        //         // pv stuff
        //         viewer.hide('*');
        //         var mol = io.pdb($scope.mutation.mutation.pdb);
        //         viewer.ballsAndSticks('mutation', mol, { color: pv.color.uniform(design.highlightColor) });
        //         viewer.ballsAndSticks('original', io.pdb($scope.mutation.original.pdb), { color: pv.color.uniform(design.lighterColor) });
        //         viewer.centerOn(mol);
        //         viewer.fitTo(mol);
        //
        //         // build graph and d3 it
        //         var residueCount = 0;
        //         var hasPreviousResidueInChain = false;
        //         var graph = { nodes : [], links : [] };
        //         $scope.mutation.original.groups.forEach(function(group) {
        //             graph.nodes.push({ "name" : group.aa + "-" + group.resn, "group" : group.aa });
        //             // connect sequential neighbors
        //             if(hasPreviousResidueInChain) {
        //                 graph.links.push({ "source" : residueCount, "target" : residueCount - 1, "value" : 3});
        //             }
        //             residueCount++;
        //             hasPreviousResidueInChain = true;
        //         });
        //
        //         // add spatially connected residue links
        //         $scope.mutation.original.hydrogenBonds.forEach(function(rri) {
        //             console.log(rri.partner1 + " -> " + rri.partner2);
        //             if(rri.partner1 > 100 || rri.partner2 > 100)
        //                 return;
        //             graph.links.push({ "source" : rri.partner1, "target" : rri.partner2, "value" : 1 });
        //         });
        //
        //         var width = 500, height = 500, radius = 2.5;
        //
        //         var color = d3.scale.category20();
        //         var force = d3.layout.force()
        //             .gravity(.05)
        //             .charge(-5)
        //             .linkDistance(3.5)
        //             .size([width, height]);
        //
        //         var svg = d3.select("#interaction-graph-original")
        //             .append("svg")
        //             .attr("width", width)
        //             .attr("height", height);
        //
        //         force.nodes(graph.nodes)
        //             .links(graph.links)
        //             .start();
        //
        //         var link = svg.selectAll(".link")
        //             .data(graph.links)
        //             .enter().append("line")
        //             .attr("class", "link")
        //             .style("stroke-width", function(d) { return Math.sqrt(d.value); });
        //
        //         var node = svg.selectAll(".node")
        //             .data(graph.nodes)
        //             .enter().append("circle")
        //             .attr("r", radius - .75)
        //             .style("fill", function(d) { return color(d.group); })
        //             .style("stroke", function(d) { return d3.rgb(color(d.group)).darker(); })
        //             .call(force.drag);
        //
        //
        //         node.append("title")
        //             .text(function(d) { return d.name; });
        //
        //         force
        //             .nodes(graph.nodes)
        //             .links(graph.links)
        //             .on("tick", tick)
        //             .start();
        //
        //         function tick() {
        //             node.attr("cx", function(d) { return d.x = Math.max(radius, Math.min(width - radius, d.x)); })
        //                 .attr("cy", function(d) { return d.y = Math.max(radius, Math.min(height - radius, d.y)); });
        //
        //             link.attr("x1", function(d) { return d.source.x; })
        //                 .attr("y1", function(d) { return d.source.y; })
        //                 .attr("x2", function(d) { return d.target.x; })
        //                 .attr("y2", function(d) { return d.target.y; });
        //         }
        //     }).catch(function(response) {
        //         console.log(response);
        //     }).finally(function () {
        //         $scope.mutationStatus = 'donso!';
        //     });
        // };

        // loading flags
        $scope.loadingModel = true;
        $scope.loadingAlignment = true;
        $scope.loadingStructure = false;

        // the full-fledged model instance
        $scope.reference = {};

        // map of pdb structures: id -> pdb-record
        $scope.chains = {};

        // the sequence alignment and variant positions therein
        $scope.alignment = {};
        var focusPosition = -1000;
        var hoverRange = [ -1000, -1000 ];
        $scope.selectionMode = false;
        $scope.windowSize = 2;

        // instance and control over pv
        $scope.features = features;
        $scope.interactionTypes = interactionTypes;
        var previousType;
        $scope.viewerOptions = {
            renderLigands : false,
            coloringFeature : features[0],
            interactionType : interactionTypes[0]
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

            // registerInViewer(chainId, true);
            initializeAlignment();
            proteinViewer = new ProteinViewerComponent('#protein-viewer');
            featureViewer = new FeatureViewerComponent('#feature-viewer', $scope.reference);
            new D3Component('#topology-viewer').topologyGraph($scope.reference);

            proteinViewer.loadStructure($scope.reference);

            console.log($scope.chains);
        }).catch(function (response) {
            console.log('model: impl error handling!');
            console.log(response);
        }).finally(function () {
            $scope.loadingModel = false;
        });

        // function registerInViewer(chainId, center) {
        //     var chain = $scope.chains[chainId];
        //     chain.pv_structure = io.pdb(chain.pdb);
        //     mol.assignHelixSheet(chain.pv_structure);
        //
        //     chain.pv_protein = viewer.cartoon('protein_' + chainId, chain.pv_structure.select('protein'), { color: pv.color.uniform(design.highlightColor) });
        //     chain.pv_ligands = viewer.ballsAndSticks('ligands_' + chainId, chain.pv_structure.select('ligand'), { color: pv.color.uniform(design.lighterColor) });
        //     if(!$scope.viewerOptions.renderLigands || $scope.selectionMode) {
        //         viewer.hide('ligands_' + chainId);
        //     }
        //     if($scope.selectionMode) {
        //         viewer.hide('protein_' + chainId);
        //         viewer.show('selection_' + chainId);
        //         // new structure registered while in selection mode: assign pv_selection property
        //         chain.pv_selection = viewer.ballsAndSticks('selection_' + chain.id,
        //             chain.pv_protein.select({ rnumRange : hoverRange }),
        //             { color: pv.color.uniform(design.highlightColor) });
        //     }
        //
        //     if(center) {
        //         viewer.centerOn(chain.pv_structure);
        //         viewer.fitTo(chain.pv_structure);
        //     }
        // }
        //
        // // fetch additional model instance from the server
        // function initializeStructure(chainId) {
        //     $scope.loadingStructure = true;
        //     // load pdb pv_structure
        //     console.log('fetching additional model for ' + chainId);
        //     ProteinService.loadModel(chainId).then(function(response) {
        //         $scope.chains[chainId] = response.data;
        //         var chain = $scope.chains[chainId];
        //         chain.selected = true;
        //         registerInViewer(chainId);
        //         // duplicated as newly fetched chains are always highlighted
        //         chain.pv_protein.setSelection(chain.pv_protein.select({ cname : chain.chainId }));
        //         chain.pv_ligands.setSelection(chain.pv_ligands.select({ cname : chain.chainId }));
        //     }).catch(function(response) {
        //         console.log('pdb: impl error handling!');
        //         console.log(response);
        //     }).finally(function() {
        //         $scope.loadingStructure = false;
        //     });
        // }
        //
        function initializeAlignment() {
            // load associated sequences from server
            ProteinService.loadAlignment($scope.reference.rep).then(function(response) {
                $scope.alignment = response.data;
                console.log($scope.alignment);

                $scope.alignment.chains.forEach(function(sequence) {
                    if(sequence.id === $scope.reference.id) {
                        $scope.reference.ec = sequence.ec;
                        $scope.reference.pfam = sequence.pfam;
                        $scope.reference.uniprot = sequence.uniprot;
                        $scope.reference.title = sequence.title;
                    }
                });

                $scope.alignment.positions.forEach(function(position) {
                    position.tooltip = position.mutant || position.variant || position.activeSite;
                });
            }).catch(function(response) {
                console.log('impl error handling at multi-sequence view:');
                console.log(response);
            }).finally(function() {
                $scope.loadingAlignment = false;
            });
        }
        //
        // /* pv controls */
        // $scope.$watch('viewerOptions.renderLigands', function(newVal, oldVal) {
        //     if(newVal === oldVal)
        //         return;
        //
        //     forAllChains(function(chain) {
        //         var id = 'ligands_' + chain.id;
        //         if(newVal && chain.selected) {
        //             viewer.show(id);
        //         } else {
        //             viewer.hide(id);
        //         }
        //     });
        //     viewer.requestRedraw();
        // });
        //
        // $scope.$watch('viewerOptions.interactionType', function(type) {
        //     if(type === interactionTypes[0])
        //         return;
        //
        //     changeVisualizedInteractions(type.tag);
        // });
        //
        // function changeVisualizedInteractions(typeTag) {
        //     if(previousType) {
        //         viewer.hide(previousType);
        //     }
        //     forAllChains(function(chain) {
        //         if (initializedInteractions.indexOf(chain + typeTag) < 0)
        //             createInteractions(chain, typeTag);
        //     });
        //     viewer.show(typeTag);
        //     viewer.requestRedraw();
        //     previousType = typeTag;
        // }
        //
        // /* render PLIP interactions */
        // function createInteractions(chain, typeTag) {
        //     console.log(chain);
        //     chain[typeTag].forEach(function(entry) {
        //         viewer.customMesh(typeTag).addTube(entry.coords1, entry.coords2, 0.15, { cap : true, color : '#FF0000' });
        //     });
        //     initializedInteractions.push(chain + typeTag);
        // }
        //
        /* interactivity functions from the msa-view */
        // $scope.enterSelectionMode = function(index) {
        //     $scope.selectionMode = true;
        //     // we update potentially once again, to ensure correct ranges in all cases
        //     focusPosition = index;
        //     determineHoverRange(index);
        //     $rootScope.context = ['selection mode', 'range: ' + hoverRange];
        //     viewer.hide('*');
        //
        //     forAllChains(function(chain) {
        //         viewer.rm('selection_' + chain.id);
        //         chain.pv_selection = viewer.ballsAndSticks('selection_' + chain.id,
        //             chain.pv_protein.select({ rnumRange : hoverRange }),
        //             { color: pv.color.uniform(design.highlightColor) });
        //         if(!chain.selected) {
        //             viewer.hide('selection_' + chain.id);
        //         }
        //         if($scope.viewerOptions.renderLigands) {
        //             viewer.show('ligands_' + chain.id);
        //         }
        //     });
        //
        //     var center = $scope.reference.pv_selection.select({ cname : $scope.reference.chainId });
        //     viewer.centerOn(center);
        //     viewer.fitTo(center);
        //     // some safety net as sometimes the highlight effect lingers after this call returns
        //     $scope.dehoverChain();
        //     viewer.requestRedraw();
        // };
        //
        // $scope.leaveSelectionMode = function() {
        //     $scope.selectionMode = false;
        //     $rootScope.context = [];
        //     forAllChains(function(chain) {
        //         viewer.rm('selection_' + chain.id);
        //         if(!chain.selected) {
        //             return;
        //         }
        //         viewer.show('protein_' + chain.id);
        //         if($scope.viewerOptions.renderLigands) {
        //             viewer.show('ligands_' + chain.id);
        //         }
        //     });
        //     viewer.centerOn($scope.reference.pv_structure.select({ cname: $scope.reference.chainId }));
        //     viewer.autoZoom();
        //     // some safety net as sometimes the highlight effect lingers after this call returns
        //     $scope.dehoverChain();
        //     viewer.requestRedraw();
        // };
        //
        // /* hover selection of residues, to be rendered as balls-and-sticks */
        // $scope.hoverPosition = function(index) {
        //     // do not visualize anything in selection mode
        //     if($scope.selectionMode) {
        //         return;
        //     }
        //
        //     focusPosition = index;
        //     determineHoverRange(index);
        //
        //     forAllChains(function(chain) {
        //         if(chain.selected === false) {
        //             return;
        //         }
        //
        //         viewer.ballsAndSticks('hover',
        //             chain.pv_protein.select({ rnumRange : hoverRange }),
        //             { color: pv.color.uniform(design.lighterColor) }
        //         );
        //     });
        // };
        //
        // function determineHoverRange(index) {
        //     hoverRange[0] = index - $scope.windowSize;
        //     if(hoverRange[0] < 0) {
        //         hoverRange[0] = 0;
        //     }
        //     hoverRange[1] = index + $scope.windowSize;
        //     if(hoverRange[1] > $scope.alignment.length) {
        //         hoverRange[1] = $scope.alignment.length;
        //     }
        // }
        //
        // $scope.isFocusPosition = function(index) {
        //     return focusPosition === index;
        // };
        //
        // $scope.isHoverPosition = function(index) {
        //     return index >= hoverRange[0] && index <= hoverRange[1];
        // };
        //
        // $scope.dehoverPosition = function() {
        //     // ignore when in selection mode
        //     if($scope.selectionMode) {
        //         return;
        //     }
        //     focusPosition = -1;
        //     hoverRange = [-1000, -1000];
        //     viewer.rm('hover');
        //     viewer.requestRedraw();
        // };
        //
        // // toggle chain in pv, potentially fetch it first
        // $scope.toggleChain = function(index) {
        //     var chainId = $scope.alignment.chains[index].id;
        //     var chain = $scope.chains[chainId];
        //     // ensure chain is loaded
        //     if(!chain || !chain.hasOwnProperty('pdb')) {
        //         initializeStructure(chainId);
        //         return;
        //     }
        //     chain.selected = !chain.selected;
        //
        //     if(chain.selected) {
        //         // either show full protein or selected fragment
        //         viewer.show(($scope.selectionMode ? 'selection_' : 'protein_') + chainId);
        //         if($scope.viewerOptions.renderLigands) {
        //             viewer.show('ligands_' + chainId);
        //         }
        //     } else {
        //         viewer.hide('protein_' + chainId);
        //         viewer.hide('ligands_' + chainId);
        //         viewer.hide('selection_' + chainId);
        //     }
        //     viewer.requestRedraw();
        // };
        //
        // $scope.hoverChain = function(index) {
        //     var chainId = $scope.alignment.chains[index].id;
        //     var chain = $scope.chains[chainId];
        //     if(!chain || !chain.selected) {
        //         return;
        //     }
        //     if($scope.selectionMode) {
        //         chain.pv_selection.setSelection(chain.pv_selection.select({ cname : chain.chainId }));
        //     } else {
        //         chain.pv_protein.setSelection(chain.pv_protein.select({ cname : chain.chainId }));
        //         chain.pv_ligands.setSelection(chain.pv_ligands.select({ cname : chain.chainId }));
        //     }
        //     viewer.requestRedraw();
        // };
        //
        // $scope.dehoverChain = function() {
        //     forAllChains(function(chain) {
        //         if(!chain.selected) {
        //             return;
        //         }
        //         var ev = chain.pv_structure.full().createEmptyView();
        //         chain.pv_protein.setSelection(ev);
        //         chain.pv_ligands.setSelection(ev);
        //         if(chain.pv_selection !== undefined) {
        //             chain.pv_selection.setSelection(ev);
        //         }
        //     });
        //     viewer.requestRedraw();
        // };
        //
        // // convenience for-all-chains function
        // function forAllChains(func) {
        //     for(var chainEntry in $scope.chains) {
        //         if(!$scope.chains.hasOwnProperty(chainEntry)) {
        //             return;
        //         }
        //         func($scope.chains[chainEntry]);
        //     }
        // }
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
            },
            mutate : function(chainId, pos, aa) {
                return $http.get(apiUrl + 'mutate/' + chainId + '/' + pos + '/' + aa);
            }
        }
    }]);

    /**
     * Formats EC numbers as link to the IUBMB.
     */
    MODULE.filter('ec', function() {
        return function(input) {
            if(input) {
                return input.split(".").join("/");
            }
        }
    });

    MODULE.filter('splitToPdbId', function() {
        return function(input) {
            if(input) {
                return input.split("_")[0];
            }
        }
    });

    MODULE.filter('mapToTlc', ['aminoAcids', function(aminoAcids) {
        return function(input) {
            if(input) {
                var tlc = input;
                aminoAcids.forEach(function(go) {
                    if(go.olc === input)
                        tlc = go.tlc;
                });
                return tlc.toUpperCase();
            }
        }
    }]);
})();
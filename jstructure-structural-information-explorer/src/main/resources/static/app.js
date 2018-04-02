'use strict';

(function () {
    var MODULE = angular.module('mpe', ['ngRoute', 'ngResource', 'ngDropdowns']);

    MODULE.constant('design', {
        defaultColor : '#454b52',
        lighterColor : '#d1e2e6',
        highlightColor : '#52b1e9'
    });

    MODULE.constant('features', [
        { name : 'average RMSD increase', tag : 'averageRmsdIncrease' },
        { name : 'average TM-score increase', tag : 'averageTmScoreIncrease' },
        { name : 'average Q increase', tag : 'averageQIncrease' },
        { name : 'maximum RMSD increase', tag : 'maximumRmsdIncrease' },
        { name : 'maximum TM-score increase', tag : 'maximumTmScoreIncrease' },
        { name : 'maximum Q increase', tag : 'maximumQIncrease' }
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

            $routeProvider.when('/chain/:stfid', {
                templateUrl: '/chain.htm',
                controller: 'ChainController'
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
    MODULE.controller('HomeController', ['$scope', '$rootScope', '$location', 'ProteinService',
        function($scope, $rootScope, $location, ProteinService) {
            $scope.ids = ProteinService.ids;
            $rootScope.prev = null;
            $rootScope.next = null;
        }]);

    /**
     * the controller for the chain-view
     */
    MODULE.controller('ChainController', ['$scope', '$rootScope', '$routeParams', 'ProteinService', 'design', 'features',
        function($scope, $rootScope, $routeParams, ProteinService, design, features) {
            var proteinViewer, featureViewer;

            // loading flags
            $scope.loadingModel = true;
            $scope.loadingStructure = false;

            // sorting of contacts/residues
            $scope.sortTypeContacts = 'residueIdentifier1';
            $scope.sortReverseContacts = false;
            $scope.sortTypeResidues = 'residueIdentifier';
            $scope.sortReverseResidues = false;

            // control over the protein-view navigation tab
            $scope.tab = 0;
            $scope.selectTab = function (setTab) {
                $scope.tab = setTab;
            };
            $scope.isSelected = function (checkTab) {
                return $scope.tab === checkTab;
            };

            $scope.featureIndex = 0;
            $scope.feature = features[$scope.featureIndex].tag;
            $scope.features = features;
            $scope.selectFeature = function (setFeatureIndex) {
                $scope.featureIndex = setFeatureIndex;
                $scope.feature = features[$scope.featureIndex].tag;
                proteinViewer.colorStructure($scope.protein, $scope.feature);
            };
            $scope.isSelectedFeature = function (checkFeatureIndex) {
                return $scope.featureIndex === checkFeatureIndex;
            };

            // init model for reference chain
            console.log('fetching model for ' + $routeParams.stfid);
            ProteinService.loadModel($routeParams.stfid).then(function (response) {
                $scope.protein = response.data;
                console.log($scope.protein);

                proteinViewer = new ProteinViewerComponent('#protein-viewer');
                featureViewer = new FeatureViewerComponent('#feature-viewer',
                    $scope.protein,
                    function(pos) {
                        proteinViewer.select(+pos.substring(0, pos.length - 1));
                    });

                proteinViewer.loadStructure($scope.protein);
                proteinViewer.colorStructure($scope.protein, $scope.feature);

                $rootScope.prev = determineLinks($scope.protein.stfId, ProteinService.ids, -1);
                $rootScope.next = determineLinks($scope.protein.stfId, ProteinService.ids, 1);
            }).catch(function (response) {
                console.log('model: impl error handling!');
                console.log(response);
            }).finally(function () {
                $scope.loadingModel = false;
            });

            function determineLinks(stfId, ids, offset) {
                var index = ids.indexOf(stfId);
                var entry = ids[index + offset];
                if(entry) {
                    return entry;
                } else {
                    return null;
                }
            }

            $scope.select = function(res1, res2) {
                proteinViewer.select(res1, res2);
            };

            $scope.deselect = function(res1, res2) {
                proteinViewer.deselect(res1, res2);
            }
        }]);

    /**
     * access to the back-end
     */
    MODULE.factory('ProteinService', ['$rootScope', '$http', function($rootScope, $http) {
        var apiUrl = '/api/';
        var ids = [];

        function logError(context, error) {
            $rootScope.alerts.push({ type: 'danger',
                msg: 'loading ' + context + ' from server failed with [' + error.status + '] ' + error.statusText
            });
        }

        $http.get(apiUrl + 'ids').then(function(response) {
            response.data.forEach(function(id) {
                ids.push(id);
            });
        }, function(response) {
            logError('all chain ids', response)
        });

        return {
            ids : ids,
            loadModel : function(chainId) {
                return $http.get(apiUrl + 'json/' + chainId);
            }
        }
    }]);
})();
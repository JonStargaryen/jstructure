'use strict';

(function () {
    var MODULE = angular.module('mpe', ['ngRoute', 'ngResource', 'nvd3', 'ngDropdowns']);

    MODULE.constant('design', {
        defaultColor : '#454b52',
        lighterColor : '#d1e2e6'
    });

    MODULE.constant('features', [
        { name : 'no coloring', tag : 'none' },
        { name : 'relative accessible surface area', tag : 'rasa' },
        { name : 'secondary structure element', tag: 'sse' }
    ]);

    MODULE.constant('renderModes',  [
        'ballsAndSticks', 'cartoon', 'sline', 'lines', 'trace', 'lineTrace', 'tube', 'spheres'
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
    MODULE.run(['$rootScope', '$http', '$location', 'renderModes', 'features', '$filter',
        function ($rootScope, $http, $location, renderModes, features, $filter) {
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

        $rootScope.features = [];
        features.forEach(function (entry) {
            $rootScope.features.push({ text: entry.name,
                index: $rootScope.features.length, raw: entry.tag });
        });
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
        $scope.protein = {};
        $scope.options = {
            renderMembrane : false,
            renderHydrogenBonds : false,
            renderMode : $scope.renderModes[0],
            coloringFeature : $scope.features[0]
        };
        var membraneInitialized = false;
        var initializedInteractions = [];

        $http.get('/api/proteins/pdbid/' + $routeParams.pdbid).then(function(data) {
            $scope.protein = data.data;

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
            visualizeProteinStructure($scope.protein);
            $scope.loading = false;
        }, function (d) {
            $rootScope.alerts.push({ type: 'danger', msg: 'loading protein ' + $routeParams.pdbid +
            ' from server failed with [' + d.status + '] ' + d.statusText });
        });

        /* PV related functions */
        var visualizeProteinStructure = function (protein) {
            var options = {
                // div to be selected is not visible at that time
                width: 400,
                height: 400,
                antialias: true,
                quality : 'high',
                background:'#313a41'
            };

            // ensure container is empty
            document.getElementById('protein-visualizer').innerHTML = '';

            $scope.viewer = pv.Viewer(document.getElementById('protein-visualizer'), options);
            $scope.structure = io.pdb($scope.protein.pdb);
            mol.assignHelixSheet($scope.structure);
            $scope.viewer.options('fog', false);

            renderStructure();

            $scope.viewer.centerOn($scope.structure);
            $scope.viewer.autoZoom();
        };

        /* render plain PDB structure with chosen style */
        var renderStructure = function() {
            $scope.viewer.renderAs('protein', $scope.structure, $scope.options.renderMode.raw,
                { color: pv.color.uniform(design.lighterColor) });
        };

        /* render ANVIL membrane layer */
        var createMembrane = function() {
            $scope.protein.membrane.forEach(function(entry) {
                $scope.viewer.customMesh('membrane').addSphere([entry[0], entry[1], entry[2]], 0.75,
                    { color : [1, 1, 1, 0.5] });
            });
            membraneInitialized = true;
        };

        /* render PLIP interactions */
        var createInteractions = function (type) {
            console.log('creating interactions for ' + type);
            $scope.protein[type].forEach(function(entry) {
                var atom1 = $scope.structure.atom(entry.a1);
                var atom2 = $scope.structure.atom(entry.a2);
                var g = $scope.viewer.customMesh(type);
                g.addTube(atom1.pos(), atom2.pos(), 0.3, { cap : false, color : interactions[type] });
            });
            initializedInteractions.push(type);
        };

        /* watchers for front-end options */
        $scope.$watch('options.renderMode', function(newVal, oldVal){
            if(newVal === oldVal)
                return;
            console.log('render mode changed from ' + oldVal + ' to ' + newVal);
            $scope.viewer.rm('protein');
            renderStructure();
        });

        $scope.$watch('options.renderMembrane', function(newVal, oldVal){
            if(newVal === oldVal)
                return;
            console.log('render membrane changed from ' + oldVal + ' to ' + newVal);
            if(newVal) {
                if(!membraneInitialized)
                    createMembrane();
                $scope.viewer.show('membrane');
            } else {
                $scope.viewer.hide('membrane');
            }
            $scope.viewer.requestRedraw();
        });

        $scope.$watch('options.renderHalogenBonds', function(newVal, oldVal){
            toggleInteractions(newVal, oldVal, 'halogenBonds');
        });

        $scope.$watch('options.renderHydrogenBonds', function(newVal, oldVal){
            toggleInteractions(newVal, oldVal, 'hydrogenBonds');
        });

        $scope.$watch('options.renderHydrophobicInteractions', function(newVal, oldVal){
            toggleInteractions(newVal, oldVal, 'hydrophobicInteractions');
        });

        $scope.$watch('options.renderMetalComplexes', function(newVal, oldVal){
            toggleInteractions(newVal, oldVal, 'metalComplexes');
        });

        $scope.$watch('options.renderPiCationInteractions', function(newVal, oldVal){
            toggleInteractions(newVal, oldVal, 'piCationInteractions');
        });

        $scope.$watch('options.renderPiStackings', function(newVal, oldVal){
            toggleInteractions(newVal, oldVal, 'piStackings');
        });

        $scope.$watch('options.renderSaltBridges', function(newVal, oldVal){
            toggleInteractions(newVal, oldVal, 'saltBridges');
        });

        $scope.$watch('options.renderWaterBridges', function(newVal, oldVal){
            toggleInteractions(newVal, oldVal, 'waterBridges');
        });

        var toggleInteractions = function(newVal, oldVal, type) {
            if(newVal === oldVal)
                return;
            console.log('render ' + type + ' changed from ' + oldVal + ' to ' + newVal);
            if(newVal) {
                if(initializedInteractions.indexOf(type) < 0)
                    createInteractions(type);
                $scope.viewer.show(type);
            } else {
                $scope.viewer.hide(type);
            }
            $scope.viewer.requestRedraw();
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
})();
'use strict';

(function () {
    var MODULE = angular.module('mpe', ['ngRoute', 'ngResource', 'nvd3', 'ngDropdowns']);

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
    MODULE.run(['$rootScope', '$http', '$location', function ($rootScope, $http, $location) {
        $rootScope.alerts = [];

        // message/alert container
        $rootScope.closeAlert = function (index) {
            $rootScope.alerts.splice(index, 1);
        };

        $rootScope.page = function () {
            return $location.path();
        }
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

    MODULE.controller('ProteinController', ['$scope', '$routeParams', '$http', '$timeout', 'ProteinService', function($scope, $routeParams, $http, $timeout, ProteinService) {
        $scope.protein = {};

        $http.get('/api/proteins/pdbid/' + $routeParams.pdbid).then(function(data) {
            $scope.protein = data.data;
        }, function (d) {
            $rootScope.alerts.push({ type: 'danger', msg: 'loading protein ' + pdbid + ' from server failed with ['+ d.status + '] '+ d.statusText });
        });
    }]);

    /**
     * fetch all available data for the home view
     */
    MODULE.controller('HomeController', ['$scope', '$location', 'ProteinService', function($scope, $location, ProteinService) {
        $scope.allProteins = ProteinService.allProteins;
        $scope.nonRedundantAlphaHelicalProteins = ProteinService.nonRedundantAlphaHelicalProteins;
    }]);

    MODULE.factory('ProteinService', ['$rootScope', '$http','$q', function($rootScope, $http, $q) {
        /* fetch all known protein ids */
        var allProteins = [];
        var nonRedundantAlphaHelicalProteins = [];

        $http.get('/api/proteins/all').then(function(d) {
            d.data.forEach(function (p) {
                allProteins.push(p);
            });
        }, function (d) {
            $rootScope.alerts.push({ type: 'danger', msg: 'loading all proteins from server failed with ['+ d.status + '] '+ d.statusText });
        });

        $http.get('/api/proteins/alpha_nr').then(function(d) {
            d.data.forEach(function (p) {
                nonRedundantAlphaHelicalProteins.push(p);
            });
        }, function (d) {
            $rootScope.alerts.push({ type: 'danger', msg: 'loading non-redundant alpha helical proteins from server failed with ['+ d.status + '] '+ d.statusText });
        });

        return { "allProteins" : allProteins,
            "nonRedundantAlphaHelicalProteins" : nonRedundantAlphaHelicalProteins };
    }]);
})();
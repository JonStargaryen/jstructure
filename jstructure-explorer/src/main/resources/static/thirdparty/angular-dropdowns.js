/**
 * @license MIT http://jseppi.mit-license.org/license.html
 */
(function(window, angular) {
    'use strict';

    var dd = angular.module('ngDropdowns', []);

    dd.run(['$templateCache', function($templateCache) {
        $templateCache.put('ngDropdowns/templates/dropdownSelect.html',
            ['<div class="wrap-dd-select" tabindex="0">',
                '<span class="selected">{{dropdownModel[labelField]}}</span>',
                '<ul class="dropdown">',
                '<li ng-repeat="item in dropdownSelect"',
                ' class="dropdown-item"',
                ' dropdown-select-item="item"',
                ' dropdown-item-label="labelField">',
                '</li>','</ul>','</div>'].join(''));

        $templateCache.put('ngDropdowns/templates/dropdownSelectItem.html',
            ['<li ',
                'ng-click="selectItem()">',
                '<a href=""',
                ' ng-href="{{dropdownMenuItem.href}}"',
                ' class="dropdown-item">',
                '{{dropdownSelectItem[dropdownItemLabel]}}',
                '</a>','</li>'].join(''));
    }]);

    dd.directive('dropdownSelect', ['DropdownService', function(DropdownService) {
        return {
            restrict : 'A',
            replace : true,
            scope : {
                dropdownSelect : '=',
                dropdownModel : '=',
                dropdownItemLabel : '@',
                dropdownOnchange : '&',
                dropdownDisabled : '='
            },

            controller : [
                '$scope',
                '$element',
                function($scope, $element) {
                    $scope.labelField = $scope.dropdownItemLabel || 'text';

                    DropdownService.register($element);

                    this.select = function(selected) {
                        if (!angular.equals(selected, $scope.dropdownModel)) {
                            $scope.dropdownModel = selected;
                        }
                        $scope.dropdownOnchange({
                            selected : selected
                        });
                        $element[0].blur(); //trigger blur to clear active
                    };

                    $element.bind('click', function(event) {
                        event.stopPropagation();
                        if (!$scope.dropdownDisabled) {
                            DropdownService.toggleActive($element);
                        }
                    });

                    $scope.$on('$destroy', function() {
                        DropdownService.unregister($element);
                    });
                }],
            templateUrl : 'ngDropdowns/templates/dropdownSelect.html'
        };
    }]);

    dd.directive('dropdownSelectItem', [ function() {
        return {
            require : '^dropdownSelect',
            replace : true,
            scope : {
                dropdownItemLabel : '=',
                dropdownSelectItem : '='
            },

            link : function(scope, element, attrs, dropdownSelectCtrl) {
                scope.selectItem = function() {
                    if (scope.dropdownSelectItem.href) {
                        return;
                    }
                    dropdownSelectCtrl.select(scope.dropdownSelectItem);
                };
            },

            templateUrl : 'ngDropdowns/templates/dropdownSelectItem.html'
        };
    }]);

    dd.factory('DropdownService', [ '$document', function($document) {
        var body = $document.find('body'), service = {}, _dropdowns = [];

        body.bind('click', function() {
            angular.forEach(_dropdowns, function(el) {
                el.removeClass('active');
            });
        });

        service.register = function(ddEl) {
            _dropdowns.push(ddEl);
        };

        service.unregister = function(ddEl) {
            var index;
            index = _dropdowns.indexOf(ddEl);
            if (index > -1) {
                _dropdowns.splice(index, 1);
            }
        };

        service.toggleActive = function(ddEl) {
            angular.forEach(_dropdowns, function(el) {
                if (el !== ddEl) {
                    el.removeClass('active');
                }
            });

            ddEl.toggleClass('active');
        };

        service.clearActive = function() {
            angular.forEach(_dropdowns, function(el) {
                el.removeClass('active');
            });
        };

        service.isActive = function(ddEl) {
            return ddEl.hasClass('active');
        };

        return service;
    }]);
})(window, window.angular);
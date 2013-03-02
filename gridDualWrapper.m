classdef gridDualWrapper < handle
   % It's just a wrapper to gridDual function
   % that supports storing number of dual
   % function calls and values returned from it.
   % If optional argument 'only_counter' is set to true,
   % then object will store only function calls counter.
   % 

   properties
      unary;
      vertC;
      horC;
      only_counter = false;
      lowerBound = [];
      upperBound = [];
      time = [];
      funcCalls = 0;
      timeStart;
   end % properties
   methods
      function obj = gridDualWrapper(unary, vertC, horC, varargin)
         obj.only_counter = process_options(varargin, 'only_counter', false);
         obj.unary = unary;
         obj.vertC = vertC;
         obj.horC = horC;
         obj.timeStart = cputime;
      end
      function [dual_energy, grad, upper_energy, labels_first, labels_second] = dual(obj, lambda)
         [dual_energy, grad, upper_energy, labels_first, labels_second] = ...
                                                   gridDual(lambda, obj.unary, obj.vertC, obj.horC);
         if (~obj.only_counter)
            obj.lowerBound(end + 1) = dual_energy;
            obj.upperBound(end + 1) = upper_energy;
            obj.time(end + 1) = cputime - obj.timeStart;
         end
         obj.funcCalls = obj.funcCalls + 1;
      end
      function [funcCalls, upperBound, lowerBound, time] = getState(obj)
         funcCalls = obj.funcCalls;
         if nargout > 1
            if ~obj.only_counter
               upperBound = obj.upperBound(:);
               lowerBound = obj.lowerBound(:);
               time = obj.time(:);
            else
               error('Object have only function calls count information!')
            end
         end
      end
   end% methods
end% classdef
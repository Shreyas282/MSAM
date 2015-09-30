function  [ Y , dtY , ddtY ] = lowessDifferentiation( t , y , varargin )
%% LOWESSDIFFERENTIATION smoothing and differentiation of scatterplots
%   [ Y ] = lowessDifferentiation( t , y ) performs smoothing on the
%           data y. The vector t describes the time of each data row in
%           y. Thus t and y must have equal number of rows, but the data 
%           y can have any number of columns. The function supports
%           non-uniform time intervals.
%
%   [ Y , dtY ] = lowessDifferentiation( t , y ) additionally output the
%           first time derivative at each data point.
%
%   [ Y , dtY , ddtY ] = lowessDifferentiation( t , y ) additionally 
%           outputs the second time derivative at each data point.
%
%   [ Y ] = lowessDifferentiation( t , y , varargin) allows for specifying
%           various parameters and options. 
%
%           The possible arguments are:
%
%               Argument        |       Value
%           ----------------------------------------------------
%           'fitpdeg'           |   any integer (default: 2)
%                               |
%           'windowhalfsize'    |   any integer (default: 7)
%                               |
%           'infnan'            |   logical (default: false)
%                               |
%           'kernel'            |   anonymous function 
%                               |       (default: @(x) ( 1 - abs(x).^3).^3)
%
%           'fitpdeg' specifies the degree of the polynomial to use for
%           fitting the data. It has a minimal value that equals the number
%           of output arguments minus one. The best value is the degree
%           that best resembles the data locally over the full window.
%
%           'windowhalfsize' specifies the smoothing interval. The optimal
%           interval depends on data under consideration. Windowhalfsize
%           has to be larger than fitpdeg/2, which is automaticcally
%           enforced and overrules the user windowhalfsize setting.
%
%           'infnan' specifies the handling of Inf and NaN values. By
%           default they are skipped but kept in place in the output. If
%           set to true, they are removed from the data. Note that in this
%           case the size of the output doesn't match the size of the input
%           and the time t doesn't match with the output for each row.
%
%           'kernel' specifies the weighting of each point in the window.
%           By default it uses tri-cube weights. You can specify any
%           (unnormalized) kernel function on [-1 1].
%

%% ALGORITHM DETAILS -----------------------------------------------------------

%   This is a fully vectorized local regression algorithm that allows you
%   to specify any weightfunction and handles missing or infinite 
%   measurements. Local regression creates a window around each points and
%   performs a weighted fit through the points in that window using a 
%   polynomial of specified degree. The implementation here constructs a 
%   linear equation for each timepoint and then solves it to find 
%   coefficients of the fitted polynomial. 
%
%   The default settings are set for default LOWESS with a windowhalfsize
%   of 7, but it is possible to create many filters using the settings.
%   For example a local equally weighted moving average is given by setting
%   'fitpdeg' to 0 and 'kernel' to (@(x) 1). 
%
%   For more details on local regression
%       https://en.wikipedia.org/wiki/Local_regression
%
%   For a list of possible kernels
%       https://en.wikipedia.org/wiki/Kernel_(statistics)
%
%   Written by Robert Noest, 2014
%

%% MAIN CODE -------------------------------------------------------------------

% Validate user input:
[fitpdeg , windowhalfsize , kernel, infnan ] = checkInput( t , y , nargout , ...
    varargin);

% Handle Inf and NaN values
if infnan % remove Inf and NaN values
    
    badRowNDX = any(~isfinite(y),2);
    
    t(badRowNDX)    = [];
    y(badRowNDX,:)  = [];
    
    [ Y , dtY , ddtY ] = runLOWESS( t , y , nargout , fitpdeg , ...
        windowhalfsize , kernel );
    
else % skip Inf and NaN rows
    
    [ Y , dtY , ddtY ] = skip( t , y , nargout , fitpdeg , windowhalfsize , ...
        kernel );
    
end

end

%% SUPPORT FUNCTIONS -----------------------------------------------------------

% Function for validating user input
function [fitpdeg , windowhalfsize , kernel, infnan ] = checkInput( t , y , numOutputs , inputArgs)

% Default settings:
FITPDEG_DEFAULT = 2;
WINDOWHALFSIZE_DEFAULT = 7;
LOWESSKERNEL_DEFAULT = @(x) ( 1 - abs(x).^3).^3;

% Check input:
if ~isvector( t )
    error('First input needs to be a vector')
end
if ~ismatrix( y )
    error('Second input needs to be a matrix')
end
if length(t) ~= size(y,1)
    error('First and second input need to have equal number of rows')
end

% Go through optional inputs
for i = 1 : 2 : length(inputArgs)
    
    if strcmpi(inputArgs{i},'fitpdeg')
        
        if length(inputArgs{i+1}(:)) == 1 && inputArgs{i+1} == round( ...
                inputArgs{i+1})
            
            fitpdeg = inputArgs{i+1};
            
        else
            
            error('fitpdeg needs to be a single integer')
            
        end
        
    elseif strcmpi(inputArgs{i},'windowhalfsize')
        
        if length(inputArgs{i+1}(:)) == 1 && inputArgs{i+1} == round( ...
                inputArgs{i+1})
            
            windowhalfsize = inputArgs{i+1};
            
        else
            
            error('Halfwidth needs to be a single integer')
            
        end
        
    elseif strcmpi(inputArgs{i},'infnan')
        
        if islogical(inputArgs{i+1}) || isnumeric(inputArgs{i+1})
            
            infnan = inputArgs{i+1}(1);
            
        else
            
            error('InfNaN setting should be a logical; true for removing')
            
        end
        
    elseif strcmpi(inputArgs{i},'kernel')
        
        if (isa(inputArgs{i+1}, 'function_handle') && strncmp(func2str(inputArgs{i+1}), '@', 1) )
            
            kernel = inputArgs{i+1};
            
        else
            
            error('Kernel should be an anonymous function')
            
        end
        
    else
        
        error('Optional arguments incorrectly specified')
        
    end
end

% Set defaults
if ~exist('fitpdeg','var')
    
    fitpdeg = FITPDEG_DEFAULT;
    
else
    
    fitpdeg = max( fitpdeg , numOutputs - 1 );
    
end

if ~exist('windowhalfsize','var')
    
    windowhalfsize = WINDOWHALFSIZE_DEFAULT;
    
else
    
    windowhalfsize = max( windowhalfsize , fitpdeg/2 + 1 );
    
end

if ~exist('kernel','var')
    
    kernel = LOWESSKERNEL_DEFAULT;
    
end

if ~exist('infnan','var')
    
    infnan = false;
    
end

end

% Function for skipping Inf and NaN values in data
function [ Y , dtY , ddtY ] = skip( t , y , numOutputs , fitpdeg , windowhalfsize , kernel )

% Copy over all data to output (for correct size and copy of NaN/Inf values)
Y = y;
if numOutputs > 1
    dtY = y;
else
    dtY = [];
end
if numOutputs > 2
    ddtY = y;
else
    ddtY = [];
end

% Find number of data columns and good data points
[ ~ , totCol ] = size(y);
goodNDX        = isfinite(y);

% Run over all columns:
colms = 1 : totCol;
while ~isempty(colms)
    
    % Find columns that have identical position of Inf/NaN values
    cols2Remove = all( bsxfun( @eq ,  goodNDX(:,colms(1)) , ...
        goodNDX(:,colms)));
    currentY = y( goodNDX(:,colms(1)) , colms(cols2Remove));
    
    % Smooth and store    
    [ tempY , tempDtY , tempDdtY ] = runLOWESS( t(goodNDX(:,colms(1))) , ...
        currentY , numOutputs , fitpdeg , windowhalfsize , kernel );
    
    currentColNDX                      = zeros(1,totCol);
    currentColNDX(colms(cols2Remove))  = 1;
    
    currentNDX = logical( bsxfun( @times, goodNDX , currentColNDX ));
    
    Y(currentNDX)   = tempY;
    if numOutputs > 1
        dtY(currentNDX) = tempDtY;
    end
    if numOutputs > 2
        ddtY(currentNDX)= tempDdtY;
    end
    
    % Remove finished columns from list
    colms(cols2Remove) = [];
    
end

end

% Function that contains the core computation
function [ Y , dtY , ddtY ] = runLOWESS( t , y , numOutputs , fitpdeg , windowhalfsize , kernel )

% Find data input size:
[ yRows , yCols ] = size(y);

% Construct array of indices of points to use inside smoothing windows:
timeIntervalNDX = bsxfun( @plus , min( max( (( 1 : yRows)' - windowhalfsize) ...
    , 1), yRows - 2 * windowhalfsize ), 0 : (2 * windowhalfsize ) ) ;

% Construct array of time values within smoothing window taking the point
% of interest as t=0:
timeWithinInterval = bsxfun( @minus , t(timeIntervalNDX) , t((1 : yRows)') ) ;

% Scale all windows such that the time domain is [-1 1]:
timeScaling = max( abs( timeWithinInterval ) , [] , 2 ) ;
timeWithinInterval = bsxfun( @rdivide , timeWithinInterval , timeScaling ) ;

% Compute all weights:
weights = kernel( timeWithinInterval );

% Compute the LHS of linear equation:
lhs = permute( sum ( bsxfun( @times , bsxfun( @power , timeWithinInterval , ...
    reshape( bsxfun( @plus , fitpdeg : -1 : 0 , ( fitpdeg : -1 : 0 )' ) , 1 ...
    , 1 , fitpdeg + 1 , fitpdeg + 1) ), weights ) , 2 ), [3 4 1 2]);

% Create array of data values within each window:
data = reshape( cell2mat( arrayfun( @( i ) y(timeIntervalNDX, i) , 1 : yCols ...
    , 'UniformOutput' , 0 )) , yRows, 2*windowhalfsize + 1 , 1 , yCols );

% Create the RHS of the linear equation:
rhs = permute( sum( bsxfun( @times , bsxfun( @power , timeWithinInterval , ...
    reshape( fitpdeg : -1 : 0 , 1 , 1 , fitpdeg + 1)) , bsxfun( @times, ...
    weights , data ) ) , 2), [3 4 1 2]);

% Solve linear system to get coefficients:
coef = permute( reshape(cell2mat(arrayfun(@(i) lhs(:,:,i) \ rhs(:,:,i) , ...
    1:yRows , 'UniformOutput' , 0 )), fitpdeg + 1, yCols , yRows), [3 2 1]);

% Read-out smoothed values from the coefficient array:
Y = coef(:, : , fitpdeg + 1);

% If available, read-out first derivative from the coefficients array:
if numOutputs > 1
    dtY = bsxfun( @rdivide , coef( :, : , fitpdeg) , timeScaling );
else
    dtY = [];
end

% If available, read-out second derivative from the coefficient array:
if numOutputs > 2;
    ddtY = bsxfun( @rdivide , 2*coef( :, : , fitpdeg - 1) , timeScaling.^2 );
else
    ddtY = [];
end

end


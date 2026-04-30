# Plotting recipes are provided by the AIBECSRecipesBaseExt and
# AIBECSRecipesParametersExt extensions. The functions below are stubs;
# methods are added when the user loads the trigger packages
# (e.g. `using Plots` for the spatial recipes, plus `using Distributions`
# for parameter PDF recipes).

function plothorizontalslice end
function plothorizontalslice! end
function surfacemap end
function surfacemap! end
function plot∫dz end
function plot∫dz! end
const plotverticalintegral = plot∫dz
const plotverticalintegral! = plot∫dz!
function plotverticalmean end
const plotverticalaverage = plotverticalmean

function plotmeridionalslice end
function plotmeridionalslice! end
function plotzonalmean end
function plotzonalmean! end
const plotzonalaverage = plotzonalmean
const plotzonalaverage! = plotzonalmean!
function plot∫dx end
const plotzonalintegral = plot∫dx

function plotzonalslice end
function plotmeridionalmean end
const plotmeridionalaverage = plotmeridionalmean
function plot∫dy end
const plotmeridionalintegral = plot∫dy

function plot∫dxdy end
function plot∫dxdy! end
const plothorizontalintegral = plot∫dxdy
const plothorizontalintegral! = plot∫dxdy!
function plothorizontalmean end
function plothorizontalmean! end
const plothorizontalaverage = plothorizontalmean
const plothorizontalaverage! = plothorizontalmean!
function plotdepthprofile end
function plotdepthprofile! end

function plottransect end

function ratioatstation end
function ratioatstation! end
function plotstencil end
function plotstencil! end

export plothorizontalslice, plothorizontalslice!, surfacemap, surfacemap!,
    plot∫dz, plot∫dz!, plotverticalintegral, plotverticalintegral!,
    plotverticalmean, plotverticalaverage,
    plotmeridionalslice, plotmeridionalslice!,
    plotzonalmean, plotzonalmean!, plotzonalaverage, plotzonalaverage!,
    plot∫dx, plotzonalintegral,
    plotzonalslice, plotmeridionalmean, plotmeridionalaverage,
    plot∫dy, plotmeridionalintegral,
    plot∫dxdy, plot∫dxdy!, plothorizontalintegral, plothorizontalintegral!,
    plothorizontalmean, plothorizontalmean!, plothorizontalaverage, plothorizontalaverage!,
    plotdepthprofile, plotdepthprofile!,
    plottransect,
    ratioatstation, ratioatstation!,
    plotstencil, plotstencil!

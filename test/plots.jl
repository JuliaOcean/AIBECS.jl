dummy = cosd.(latvec(grd)) + sqrt.(depthvec(grd)) / 30
xs = (dummy, dummy*u"mol/m^3")

nokwplots = (surfacemap, plotverticalintegral, plotverticalaverage, plotzonalaverage, plotzonalintegral, plotmeridionalaverage, plotmeridionalintegral, plothorizontalintegral, plothorizontalaverage)
@testset "plots without kwargs" begin
    @testset "$f(x, grd)" for f in nokwplots, x in xs
        @test f(x, grd) isa Plots.Plot
    end
end

@testset "Plots with kwargs" begin
    lons = (330, 330°)
    lats = (30, 30°)
    depths = (500, 500m, 5km)
    lonlats = ((-30,30), (30,30), (-30°,30°))
    @testset "plothorizontalslice depth=$depth" for x in xs, depth in depths
        @test plothorizontalslice(x, grd, depth=depth) isa Plots.Plot
    end
    @testset "plotzonalslice lat=$lat" for x in xs, lat in lats
        @test plotzonalslice(x, grd, lat=lat) isa Plots.Plot
    end
    @testset "plotmeridionalslice lon=$lon" for x in xs, lon in lons
        @test plotmeridionalslice(x, grd, lon=lon) isa Plots.Plot
    end
    @testset "plotdepthprofile lonlat=$lonlat" for x in xs, lonlat in lonlats
        @test plotdepthprofile(x, grd, lonlat=lonlat) isa Plots.Plot
        plt = plotdepthprofile(x, grd, lonlat=lonlat)
        @test plotdepthprofile!(plt, x, grd, lonlat=lonlats[2]) isa Plots.Plot
    end
end




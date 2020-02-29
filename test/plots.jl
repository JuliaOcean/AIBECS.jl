dummy = cosd.(latvec(grid)) + sqrt.(depthvec(grid)) / 30
xs = (dummy, dummy*u"mol/m^3")

nokwplots = (surfacemap, verticalaverage, verticalintegral, verticalaverage, zonalaverage)

@testset "plots without kwargs" begin
    @testset "$f(x, grid)" for f in nokwplots, x in xs
        @test f(x,grid) isa Plots.Plot
    end
end

@testset "Plots with kwargs" begin
    lons = (330, 330u"°")
    lats = (30, 30u"°")
    depths = (500, 500u"m", 5u"km")
    depthlims = ((500,1000), (500, 1000).*u"m", (0.5,1).*u"km")
    lonlats = ((-30,30), (30, 30), (-30, 30).*u"°")
    masks = (latvec(grid) .< -40, 200 .< lonvec(grid) .< 300)

    @testset "horizontalslice depth=$depth" for x in xs, depth in depths
        @test horizontalslice(x, grid, depth=depth) isa Plots.Plot
    end
    @testset "verticalintegral depthlim=$depthlim" for x in xs, depthlim in depthlims
        @test verticalintegral(x, grid, depthlim=depthlim) isa Plots.Plot
    end
    @testset "verticalaverage depthlim=$depthlim" for x in xs, depthlim in depthlims
        @test verticalaverage(x, grid, depthlim=depthlim) isa Plots.Plot
    end
    @testset "zonalslice lat=$lat" for x in xs, lon in lons
        @test zonalslice(x, grid, lat=lat) isa Plots.Plot
    end
    @testset "meridionalslice lon=$lon" for x in xs, lat in lats
        @test meridionalslice(x, grid, lon=lon) isa Plots.Plot
    end
    @testset "zonalaverage with mask" for x in xs, mask in masks
        @test zonalaverage(x, grid, mask=mask) isa Plots.Plot
    end
    @testset "depthprofile lonlat=$lonlat" for x in xs, lonlat in lonlats
        @test depthprofile(x, grid, lonlat=lonlat) isa Plots.Plot
    end


end




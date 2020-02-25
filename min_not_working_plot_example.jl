using PyPlot


lines = [[[(10,30), (20,70)] [(20,40),(30, 100)]], [[(1,5), (2,7)] [(2,4),(3, 1)]]]


ylims = 

fig = figure()
ax = gca()


for i=1:2

    ls = matplotlib.collections.LineCollection(lines[i])
    ax.clear()
    ax.add_collection(ls)
    gui()
    #show()
    #axis("image")
    sleep(3)
end


# Having Plot as an argument is a lazy way to avoid using RecipesBase
function gray_square(xs,ys,Plot; inside = true)
  if inside == true
    dx = abs(xs[2]- xs[1])/160.0;
    dy = abs(ys[2]- ys[1])/160.0;
    ys[1] += dy
    ys[2] -= dy
    xs[1] += dx
    xs[2] -= dx
  end
    Plot([xs[1],xs[2]],[ys[2],ys[2]], linecolor = :grey, linestyle=:dash, lab="",border=false);
    Plot([xs[2],xs[2]],[ys[2],ys[1]], linecolor = :grey, linestyle=:dash, lab="",border=false);
    Plot([xs[2],xs[1]],[ys[1],ys[1]], linecolor = :grey, linestyle=:dash, lab="",border=false);
    Plot([xs[1],xs[1]],[ys[1],ys[2]], linecolor = :grey, linestyle=:dash, lab="",border=false)
end

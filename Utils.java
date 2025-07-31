package com.mygame;

import com.badlogic.gdx.Gdx;
import com.badlogic.gdx.graphics.Color;
import com.badlogic.gdx.graphics.OrthographicCamera;
import com.badlogic.gdx.graphics.glutils.ShapeRenderer;
import java.util.Iterator;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Arrays;
import java.util.List;
public class Utils<T extends Utils.Entity<T>>{

    boolean canRend = false;
    OrthographicCamera camera;
    ShapeRenderer rend;
    Quadtree<T> master;
    LinkedList<T> points = new LinkedList<T>();
    int max;
    Vec  tl = N(), br = N(Gdx.graphics.getWidth(), Gdx.graphics.getHeight());
    
    public static class Vec {
        double x;
        double y;
        private Vec(double x, double y) {
            this.x = x;
            this.y = y;
        }
    }
    
    public Utils(OrthographicCamera camera, ShapeRenderer rend, int max) {
        this.camera = camera;
        this.rend = rend;
        this.canRend = true;
        this.max = max;
        master = new Quadtree<T>(null, tl, br, max, this);
    }

    public Utils() {

    }

    public Vec N(Vec v) {
        return new Vec(v.x, v.y);
    }

    public Vec N(double x, double y) {
        return new Vec(x,y);
    }

    public Vec N(int x, int y) {
        return new Vec((double)x,(double)y);
    }

    public Vec N() {
        return new Vec(0.0,0.0);
    }
    public Vec NR() {
        return new Vec(Math.random() * 0.1, Math.random() * 0.1);
    }

    public double polygonarea(Vec[] hull) {
        double area = 0;
        for (int i = 0; i < hull.length; i++) {
            Vec p1 = hull[i];
            Vec p2 = hull[(i + 1) % hull.length];
            area += (p1.x * p2.y - p2.x * p1.y);
        }
        return area * 0.5;
    }

    public double Cross3(Vec p1, Vec p2, Vec p3) {
        return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
    }

    public Vec Sub(Vec p1, Vec p2) {
        return new Vec(p1.x - p2.x, p1.y - p2.y);
    }

    public Vec Add(Vec p1, Vec p2) {
        return new Vec(p1.x + p2.x, p1.y + p2.y);
    }
    
    public Vec Mul(Vec p, double s) {
        return new Vec(p.x * s, p.y * s);
    }

    public Vec Div(Vec p, double s) {
        return new Vec(p.x / s, p.y / s);
    }

    public Vec Permul(Vec p1, Vec p2) {
        return new Vec(p1.x * p2.x, p1.y * p2.y);
    }
    public Vec Abs(Vec v) {
        return new Vec(Math.abs(v.x), Math.abs(v.y));
    }
    public Vec Max(Vec p1, Vec p2) {
        return N(Math.max(p1.x, p2.x),Math.max(p1.y,p2.y));
    }
    public Vec Min(Vec p1, Vec p2) {
        return N(Math.min(p1.x, p2.x),Math.min(p1.y,p2.y));
    }
    public double Dot(Vec p1, Vec p2) {
        return (p1.x * p2.x +  p1.y * p2.y);
    }
    public double Area(Vec v) {
        return v.x * v.y;
    }
     public boolean LessThan(Vec p, Vec v) {
        return (p.x < v.x && p.y < v.y) ? true : false;
    }
    public boolean MoreThan(Vec p, Vec v) {
        return (p.x > v.x && p.y > v.y) ? true : false;
    }
    public double Det(Vec p1, Vec p2) {
        return (p1.x * p2.y -  p1.y * p2.x);
    }
    
    public boolean Eq(Vec p1, Vec p2) {
        return ((p1.x == p2.x) && (p1.y == p2.y));
    }

    public double Len2(Vec p1) {
        return (p1.x * p1.x + p1.y * p1.y);
    }
    public double Len(Vec p1) {
        return Math.sqrt(p1.x * p1.x + p1.y * p1.y);
    }

    public double Dist(Vec p1, Vec p2) {
        return Len(Sub(p1, p2));
    }

    public Vec Norm(Vec p1) {
        double len = Len(p1);
        return Div(p1, (len != 0) ? len : 1);
    }

    public Vec Rotate(double angle, Vec p1) {
        double x = p1.x, y = p1.y;
        double c = Math.cos(angle), s = Math.sin(angle);
        return new Vec(c*x - s*y, s*x + c*y);
    }

    public double smoothstep(double edge0, double edge1, double x) {
        double t = Clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
        return t*t * (3 - 2 * t);
    }
    public Vec Perp(Vec p1) {
        return new Vec(-p1.y, p1.x);
    }
    public Vec RevPerp(Vec p1) {
        return new Vec(p1.y, -p1.x);
    }

    public Vec Clamp(Vec p, Vec min, Vec max) {
        return new Vec(Math.max(Math.min(max.x, p.x), min.x), Math.max(Math.min(max.y, p.y), min.y));
    }

    public double Clamp(double num, double min, double max) {
        return Math.max(Math.min(max, num), min);
    }

    

//shapes
    public class ConvexHull {
        List<T> hull;
        public ConvexHull(T[] shape) {
            this.hull = update(shape);
        }

        public void DrawHull(Color color) {
            if (this.hull == null || this.hull.size() < 3) return;
            rend.setProjectionMatrix(camera.combined);
            rend.begin(ShapeRenderer.ShapeType.Line);
            rend.setColor(color);
            for (int i = 0; i < this.hull.size(); i++) {
                int j = (i == 0) ? this.hull.size() -1 : i - 1;
                T p1 = this.hull.get(i);
                T p2 = this.hull.get(j);
                rend.line((float)p1.getPos().x, (float)p1.getPos().y, (float)p2.getPos().x, (float)p2.getPos().y);
                
            }
        }

        public void drawHullNormals(ShapeRenderer renderer) {
            Vec center = this.center();
            int count = this.hull.size();

            renderer.begin(ShapeRenderer.ShapeType.Line);
            renderer.setColor(Color.RED);

            for (int i = 0; i < count; i++) {
                int j = (i == 0) ? count - 1 : i - 1;
                Vec a = this.hull.get(i).getPos();
                Vec b = this.hull.get(j).getPos();

                Vec edge = Sub(b, a);
                Vec mid = Mul(Add(a, b), 0.5);
                Vec norm = RevPerp(edge);

                Vec toCenter = Sub(center, a);
                if (Dot(norm, toCenter) > 0) {
                    norm = Mul(norm, -1);
                }

                Vec normalDir = Norm(norm);
                Vec normalEnd = Add(mid, Mul(normalDir, 10)); // 10 pixels length

                renderer.line((float)mid.x, (float)mid.y, (float)normalEnd.x, (float)normalEnd.y);
            }

            renderer.end();
        }

        private void pop(ArrayList<T> list) {
            list.remove(list.size() - 1);
        }

        public List<T> update(T[] shape) {
            T[] sortedShape = shape.clone();

            Arrays.sort(sortedShape, (a,b) -> {
                return (a.getPos().x != b.getPos().x) ? Double.compare(a.getPos().x, b.getPos().x) : Double.compare(a.getPos().y, b.getPos().y);
            });

            ArrayList<T> bottom = new ArrayList<>();
            ArrayList<T> top = new ArrayList<>();

            for (int i = 0; i < sortedShape.length; i++) {
                if (bottom.size() < 3) {
                    bottom.add(sortedShape[i]);
                } else {
                    while (bottom.size() >= 3 && 0 >= Cross3(bottom.get(bottom.size() - 2).getPos(), bottom.get(bottom.size() - 1).getPos(), sortedShape[i].getPos())) {
                        pop(bottom);
                    } 
                    bottom.add(sortedShape[i]);
                }
            }

            for (int i = sortedShape.length - 1; i >= 0 ; i--) {
                if (top.size() < 3) {
                    top.add(sortedShape[i]);
                } else {
                    while (top.size() >= 3 && 0 >= Cross3(top.get(top.size() - 2).getPos(), top.get(top.size() - 1).getPos(), sortedShape[i].getPos()) ) {
                        pop(top);
                    } 
                    top.add(sortedShape[i]);
                }
            }

            bottom.remove(bottom.size() - 1);
            top.remove(top.size() - 1);
            bottom.addAll(top);
            
            return bottom;
        }   
        public Vec center() {
            Vec v = N();
            for (T vert : this.hull) {
                v = Add(v, vert.getPos());
            }
            return Div(v, this.hull.size());
        }
        public Vec PIP(T v)  {
            Vec mtv = null;
            double MP = Double.MAX_VALUE;

            int count = this.hull.size();
 
            for (int i = 0; i < count; i++ ) {
                int j = (i == 0) ? count - 1 : i - 1;
                Vec a = this.hull.get(i).getPos();
                Vec b = this.hull.get(j).getPos();

                Vec edge = Sub(b, a);

                Vec norm = Perp(edge);

                Vec center = this.center(); // approximate shape center
                Vec toPoint = Sub(v.getPos(), center);

                if (Dot(norm, toPoint) < 0) {
                    norm = Mul(norm, -1); // flip to face outward
                }

                double minP = Double.MAX_VALUE;
                double maxP = -Double.MAX_VALUE;

                for (T vert : this.hull) {
                    double proj = Dot(vert.getPos(), norm);
                    minP = Math.min(minP, proj);
                    maxP = Math.max(maxP, proj);
                }

                double pointProj = Dot(v.getPos(), norm);
 
                if (pointProj < minP || pointProj > maxP) {
                    return null;
                }
                double pen = maxP - pointProj; // assumes point is always inside

                if (pen < MP) {
                    MP = pen;
                    mtv = norm;
                }
                
            }
            return (mtv != null) ? Mul(mtv, MP) : null;
        }
    }
    public class Circle {
        Vec center;
        double rad;
        double area;
        Circle(Vec center, double rad) {
            this.center = center;
            this.rad = rad;
            this.area = 2 * Math.PI * rad * rad;
        }
    }
    
    public class Rectangle {
        Vec start, end, bounds;
        double area;
        public Rectangle(Vec start, Vec end, Vec bounds, double area) {
            this.start = start;
            this.end = end;
            this.bounds = bounds;
            this.area = area;
        }
        public Rectangle(Vec start, Vec end) {
            this.start = start;
            this.end = end;
            this.bounds = Sub(end, start);
            this.area = Area(this.bounds);
        }

        public void MoveRec(Rectangle rec, Vec v) {
            rec.start = Add(rec.start, v);
            rec.end = Add(rec.end, v);
        }

        public void Update(Vec start, Vec end) {
            this.start = (start != null) ? start : this.start;
            this.end = (end != null) ? end : this.end;
            this.bounds = Sub(this.end, this.start);
            this.area = Area(this.bounds);
        }

        public Rectangle Copy(Rectangle rect) {
            return new Rectangle(this.start, this.end, this.bounds, this.area);
        }
    }

    public Rectangle difference(Rectangle a, Rectangle b) {
        Vec overlap = Max(a.start,b.start);
        Vec overlapB = Sub(Min(Add(a.start, a.bounds),Add(b.start, b.bounds)), overlap);
        if (overlap.x > 0 && overlap.y > 0) {
            return new Rectangle(overlap, Add(overlap, overlapB));
        } else {
            return null;
        }
    }

    public Double Area(Vec... Vecs) {
        double area = 0;
        for (int i = 0; i < Vecs.length - 1; i++ ) {
            area += Det(Vecs[i], Vecs[i + 1]);
        }
        area += Det(Vecs[Vecs.length], Vecs[0]);
        return area;
    }

    //quadtree stuff

    

    public void DrawCircs() {
        for (T i : this.points) {
            if (i.getRad() != null) i.Draw(rend);
        }
    }

    public void Insert(T v) {
        points.add(v);
        master.AddQuad(points.get(points.size()-1));
    }

    public void Refresh() {
        Iterator<T> it = points.iterator();
        ArrayList<Quadtree<T>> cull = new ArrayList<>();
        while (it.hasNext()) {
            T start = it.next();

            Quadtree<T> parent = start.getPar();
            if (parent == null) {System.out.println("orphan"); continue;} 
            
            parent.Vecs.remove(start);
            master.AddQuad(start);
            int total = 0;
            if (parent != null && parent.Parent != null) {
                List<Quadtree<T>> siblings = parent.Parent.children;
                if (siblings == null) continue;
                // do stuff
                for (Quadtree<T> child : siblings) {
                    if (child == null) continue;
                    if (child.Vecs != null) {
                        total += child.Vecs.size();
                    }
                    
                }
                if (total < this.max / 2) {
                    if (parent == null || parent.Parent == null || siblings == null) continue;
                    cull.add(parent.Parent);
                }
            }
                
            
        }
        Cull(cull);
    }

    public void Cull(List<Quadtree<T>> layers) {
        for (Quadtree<T> Q : layers) {
            if (Q.children == null) continue;
            for (Quadtree<T> child : Q.children) {
                if (child == null) continue;
                for (int i = 0; i < child.Vecs.size(); i++) {
                    Q.Vecs.add(child.Vecs.get(i));
                    child.Vecs.get(i).setPar(Q);
                }
            }
            Q.children.clear();;
            Q.isLeaf = true;
        }
    }
    
    public void ScreenResize(Vec br) {
        master = new Quadtree<T>(null, N(), br, max, this);
        for (int i  = 0; i < points.size(); i++) {
            master.AddQuad(points.get(i));
        }
    }

    public List<T> BroadSearch(Vec p1, Vec p2) {
        List<Quadtree<T>> valid = master.ValidQuadrants(p1, p2);
        ArrayList<T> found = new ArrayList<>();
        for (Quadtree<T> grid : valid) {
            for (T e : grid.Vecs) {
                found.add(e);
            }
        }
        return found;
    }

    public List<T> CircleSearch(T e) {
        List<Quadtree<T>> valid = master.ValidQuadrants(e);
        ArrayList<T> found = new ArrayList<>();
        for (Quadtree<T> grid : valid) {
            for (T p : grid.Vecs) {
                Vec d = Sub(e.getPos(), p.getPos());
                double tr = e.getRad() + p.getRad();

                if (Len2(d) <= tr * tr) {
                    found.add(p);
                }
            }
        }
        return found;
    }
    


    public interface Entity<T extends Entity<T>> {
        public Vec getPos();
        public Double getRad();
        public Utils.Quadtree<T> getPar();
        public void setPar(Utils.Quadtree<T> par);

        default void Draw(ShapeRenderer renderer) {
            Double r = getRad();
            if (r != null) {
                Vec p = getPos();
                renderer.circle((float)p.x, (float)p.y, (float)(double)r);
            }
        }
        
        default void Update() {
            System.out.println("");
        }
    }

    public static class Quadtree<T extends Entity<T>> {
        int size;
        int max;
        Quadtree<T> Parent;
        Vec b1, b2;
        boolean isLeaf = true;
        Vec dif;
        ArrayList<T> Vecs = new ArrayList<>();
        List<Quadtree<T>> children = new ArrayList<>(4);
        Utils<T> Vec2;

        public Quadtree(Quadtree<T> parent, Vec topLeft, Vec bottomRight, int max, Utils<T> u) {
            this.Vec2 = u;
            this.Parent = parent;
            this.b1 = topLeft;
            this.b2 = bottomRight;
            this.dif = this.Vec2.Sub(this.b2, this.b1);
            this.max = max;
        }

        public void AddQuad(T v) {
            if (this.isLeaf) {
                v.setPar(this);
                Vecs.add(v);
                if (Vecs.size() > this.max) {
                    List<T> old = Vecs;
                    Vecs = new ArrayList<>();
                    subdivide();
                    for (int i = 0; i < old.size(); i++) {
                        this.AddQuad(old.get(i));
                    }
                    this.isLeaf = false;
                }
            } else {
                boolean Added = false;
                for (int i = 0; i < 4 && !Added; i++) {
                    Quadtree<T> Child = this.children.get(i);

                    if (Child.Contains(v)) {
                        Child.AddQuad(v);
                        return;
                    }
                }
                v.setPar(this);
                Vecs.add(v);
            }
        }

        public boolean Contains(T v) {
            return (this.Vec2.LessThan(v.getPos(), this.b2) && (this.Vec2.MoreThan(v.getPos(), this.b1) || this.Vec2.Eq(v.getPos(), this.b1))) ? true : false;
        }

        private boolean overlaps(Vec p1, Vec p2) {
            return !(p2.x < this.b1.x || p1.x > this.b2.x || p2.y < this.b1.y || p1.y > this.b2.y);
        }

        public boolean overlaps(T e) {
            if (e.getRad() == null) {
                System.out.println("failed");
                return false;
            }
            double minX = Math.min(this.b1.x, this.b2.x);
            double maxX = Math.max(this.b1.x, this.b2.x);
            double minY = Math.min(this.b1.y, this.b2.y);
            double maxY = Math.max(this.b1.y, this.b2.y);

            double clampedX = Math.max(minX, Math.min(e.getPos().x, maxX));
            double clampedY = Math.max(minY, Math.min(e.getPos().y, maxY));
            Vec n = new Vec(clampedX, clampedY);

            Vec d = this.Vec2.Sub(e.getPos(), n);
            System.out.println(this.Vec2.Len2(d) <= e.getRad() * e.getRad());
            return this.Vec2.Len2(d) <= e.getRad() * e.getRad();
        }


        public void Draw(ShapeRenderer rend) {
            
            rend.begin(ShapeRenderer.ShapeType.Line);
            
            Vec t = this.Vec2.Sub(b2, b1);
            rend.setColor(Color.BLUE);
            rend.rect((float)this.b1.x, (float)this.b1.y,(float) t.x, (float)t.y);
            rend.end();
            if (!this.isLeaf) {
                for (int i = 0; i <  4; i++) {
                    this.children.get(i).Draw(rend);
                }
            }
        }

        public void subdivide() {
            this.isLeaf = false;
            Vec mid = this.Vec2.Div(this.Vec2.Sub(b2, b1), 2);
            Vec midPoint = this.Vec2.Add(b1, mid);

            this.children.clear();
            this.children.addAll(List.of(new Quadtree<T>(this, b1, midPoint, max, this.Vec2),

                // Top-right
                new Quadtree<T>(this, new Vec(midPoint.x, b1.y), new Vec(b2.x, midPoint.y), max,this.Vec2),

                // Bottom-left
                new Quadtree<T>(this, new Vec(b1.x, midPoint.y), new Vec(midPoint.x, b2.y), max,this.Vec2),

                // Bottom-right
                new Quadtree<T>(this, midPoint, b2, max,this.Vec2)));
        }
        public ArrayList<Quadtree<T>> ValidQuadrants(Vec p1, Vec p2) {
            Vec dif = this.Vec2.Sub(p2, p1);
            ArrayList<Quadtree<T>> valid = new ArrayList<>();
            return ValidQuadrants(p1,p2,dif,valid);
        }
        public ArrayList<Quadtree<T>> ValidQuadrants(Vec p1, Vec p2, Vec dif, ArrayList<Quadtree<T>> current) {
            if (!overlaps(p1, p2)) return current;

            if (this.isLeaf) {
                current.add(this);
            } else {
                for (Quadtree<T> child : children) {
                    child.ValidQuadrants(p1, p2, dif, current);
                }
            }
            return current;
        }
        public ArrayList<Quadtree<T>> ValidQuadrants(T e) {
            ArrayList<Quadtree<T>> valid = new ArrayList<>();
            return ValidQuadrants(e,valid);
        }

        public ArrayList<Quadtree<T>> ValidQuadrants(T e, ArrayList<Quadtree<T>> current) {
            if (this.overlaps(e) == false) { return current;}

            if (this.isLeaf) {
                current.add(this);
            } else {
                for (Quadtree<T> child : children) {
                    child.ValidQuadrants(e, current);
                }
            }
            return current;
        }
    }
}

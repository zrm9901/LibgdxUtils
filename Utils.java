package com.mygame;

import com.badlogic.gdx.graphics.Color;
import com.badlogic.gdx.graphics.OrthographicCamera;
import com.badlogic.gdx.graphics.glutils.ShapeRenderer;
import java.util.Iterator;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Arrays;
public class Utils{
    boolean canRend = false;
    OrthographicCamera camera;
    ShapeRenderer rend;
    public static class Vec {
        double x;
        double y;
        private Vec(double x, double y) {
            this.x = x;
            this.y = y;
        }

    }
    
    public Utils(OrthographicCamera camera, ShapeRenderer rend) {
        this.camera = camera;
        this.rend = rend;
        this.canRend = true;
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
        Entity[] hull;
        public ConvexHull(Entity[] shape) {
            this.hull = update(shape);
        }

        public void DrawHull(Color color) {
            if (this.hull == null || this.hull.length < 3) return;
            rend.setProjectionMatrix(camera.combined);
            rend.begin(ShapeRenderer.ShapeType.Line);
            rend.setColor(color);
            for (int i = 0; i < this.hull.length; i++) {
                int j = (i == 0) ? this.hull.length -1 : i - 1;
                Entity p1 = this.hull[i];
                Entity p2 = this.hull[j];
                rend.line((float)p1.pos.x, (float)p1.pos.y, (float)p2.pos.x, (float)p2.pos.y);
                
            }
        }

        public void drawHullNormals(ShapeRenderer renderer) {
            Vec center = this.center();
            int count = this.hull.length;

            renderer.begin(ShapeRenderer.ShapeType.Line);
            renderer.setColor(Color.RED);

            for (int i = 0; i < count; i++) {
                int j = (i == 0) ? count - 1 : i - 1;
                Vec a = this.hull[i].pos;
                Vec b = this.hull[j].pos;

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

        private void pop(ArrayList<Entity> list) {
            list.remove(list.size() - 1);
        }

        public Entity[] update(Entity[] shape) {
            Entity[] sortedShape = shape.clone();

            Arrays.sort(sortedShape, (a,b) -> {
                return (a.pos.x != b.pos.x) ? Double.compare(a.pos.x, b.pos.x) : Double.compare(a.pos.y, b.pos.y);
            });

            ArrayList<Entity> bottom = new ArrayList<>();
            ArrayList<Entity> top = new ArrayList<>();

            for (int i = 0; i < sortedShape.length; i++) {
                if (bottom.size() < 3) {
                    bottom.add(sortedShape[i]);
                } else {
                    while (bottom.size() >= 3 && 0 >= Cross3(bottom.get(bottom.size() - 2).pos, bottom.get(bottom.size() - 1).pos, sortedShape[i].pos)) {
                        pop(bottom);
                    } 
                    bottom.add(sortedShape[i]);
                }
            }

            for (int i = sortedShape.length - 1; i >= 0 ; i--) {
                if (top.size() < 3) {
                    top.add(sortedShape[i]);
                } else {
                    while (top.size() >= 3 && 0 >= Cross3(top.get(top.size() - 2).pos, top.get(top.size() - 1).pos, sortedShape[i].pos) ) {
                        pop(top);
                    } 
                    top.add(sortedShape[i]);
                }
            }

            bottom.remove(bottom.size() - 1);
            top.remove(top.size() - 1);
            bottom.addAll(top);
            
            return bottom.toArray(new Entity[0]);
        }   
        public Vec center() {
            Vec v = N();
            for (Entity vert : this.hull) {
                v = Add(v, vert.pos);
            }
            return Div(v, this.hull.length);
        }
        public Vec PIP(Entity v)  {
            Vec mtv = null;
            double MP = Double.MAX_VALUE;

            int count = this.hull.length;
 
            for (int i = 0; i < count; i++ ) {
                int j = (i == 0) ? count - 1 : i - 1;
                Vec a = this.hull[i].pos;
                Vec b = this.hull[j].pos;

                Vec edge = Sub(b, a);

                Vec norm = Perp(edge);

                Vec center = this.center(); // approximate shape center
                Vec toPoint = Sub(v.pos, center);

                if (Dot(norm, toPoint) < 0) {
                    norm = Mul(norm, -1); // flip to face outward
                }

                double minP = Double.MAX_VALUE;
                double maxP = -Double.MAX_VALUE;

                for (Entity vert : this.hull) {
                    double proj = Dot(vert.pos, norm);
                    minP = Math.min(minP, proj);
                    maxP = Math.max(maxP, proj);
                }

                double pointProj = Dot(v.pos, norm);
 
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


    public static class Entity {
        Vec pos;
        Quadtree parentNode;
        int index;
        Double rad;

        public Entity(Vec pos) {
            this.pos = pos;
        }
        public Entity(Vec pos, Quadtree parent) {
            this.pos = pos;
            this.parentNode = parent;
        }
        public Entity(Vec pos, Quadtree parent, int index) {
            this.pos = pos;
            this.parentNode = parent;
            this.index = index;
        }
        public Entity(Vec pos, double rad) {
            this.pos = pos;
            this.rad = rad;
        }
        public Entity(Vec pos, Quadtree parent, double rad) {
            this.pos = pos;
            this.parentNode = parent;
            this.rad = rad;
        }
        public Entity(Vec pos, Quadtree parent, int index, double rad) {
            this.pos = pos;
            this.parentNode = parent;
            this.index = index;
            this.rad = rad;
        }

        public void Draw(ShapeRenderer renderer) {
            if (this.rad == null) return;
            double d = this.rad;
            renderer.circle((float)this.pos.x, (float)this.pos.y, (float)d);
        }
    }

    public class QuadWorld {
        Quadtree master;
        LinkedList<Entity> points = new LinkedList<Entity>();
        int max;

        public QuadWorld(Vec topLeft, Vec bottomRight, int max) {
            this.max = max;
            master = new Quadtree(null, N(), bottomRight, max);
        }

        public void DrawCircs() {
            for (Entity i : this.points) {
                if (i.rad != null) i.Draw(rend);
            }
        }

        public void Insert(Entity v) {
            points.add(v);
            master.AddQuad(points.get(points.size()-1));
        }

        public void Refresh() {
            Iterator<Entity> it = points.iterator();
            ArrayList<Quadtree> cull = new ArrayList<>();
            while (it.hasNext()) {
                Entity start = it.next();

                Quadtree parent = start.parentNode;
                if (parent == null) {System.out.println("orphan"); continue;} 
                
                parent.Vecs.remove(start);
                master.AddQuad(start);
                int total = 0;
                if (parent != null && parent.Parent != null) {
                    Quadtree[] siblings = parent.Parent.children;
                    if (siblings == null) continue;
                    // do stuff
                    for (Quadtree child : siblings) {
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
            Cull(cull.toArray(new Quadtree[0]));
        }

        public void Cull(Quadtree[] layers) {
            for (Quadtree Q : layers) {
                if (Q.children == null) continue;
                for (Quadtree child : Q.children) {
                    if (child == null) continue;
                    for (int i = 0; i < child.Vecs.size(); i++) {
                        Q.Vecs.add(child.Vecs.get(i));
                        child.Vecs.get(i).parentNode = Q;
                    }
                }
                Q.children = new Quadtree[4];
                Q.isLeaf = true;
            }
        }
        
        public void ScreenResize(Vec br) {
            master = new Quadtree(null, N(), br, max);
            for (int i  = 0; i < points.size(); i++) {
                master.AddQuad(points.get(i));
            }
        }

        public Entity[] BroadSearch(Vec p1, Vec p2) {
            Quadtree[] valid = master.ValidQuadrants(p1, p2).toArray(new Quadtree[0]);
            ArrayList<Entity> found = new ArrayList<>();
            for (Quadtree grid : valid) {
                for (Entity e : grid.Vecs) {
                    found.add(e);
                }
            }
            return found.toArray(new Entity[0]);
        }

        public Entity[] CircleSearch(Entity e) {
            Quadtree[] valid = master.ValidQuadrants(e).toArray(new Quadtree[0]);
            ArrayList<Entity> found = new ArrayList<>();
            for (Quadtree grid : valid) {
                for (Entity p : grid.Vecs) {
                    Vec d = Sub(e.pos, p.pos);
                    double tr = e.rad + p.rad;

                    if (Len2(d) <= tr * tr) {
                        found.add(p);
                    }
                }
            }
            return found.toArray(new Entity[0]);
        }
    }

    public class Quadtree {
        int size;
        int max;
        Quadtree Parent;
        Vec b1, b2;
        boolean isLeaf = true;
        Vec dif;
        ArrayList<Entity> Vecs = new ArrayList<>();
        Quadtree[] children = new Quadtree[4];

        public Quadtree(Quadtree parent, Vec topLeft, Vec bottomRight, int max ) {
            this.Parent = parent;
            this.b1 = topLeft;
            this.b2 = bottomRight;
            this.dif = Sub(this.b2, this.b1);
            this.max = max;
        }

        public void AddQuad(Entity v) {
            if (this.isLeaf) {
                v.parentNode = this;
                Vecs.add(v);
                if (Vecs.size() > this.max) {
                    Entity[] old = Vecs.toArray(new Entity[0]);
                    Vecs = new ArrayList<>();
                    subdivide();
                    for (int i = 0; i < old.length; i++) {
                        this.AddQuad(old[i]);
                    }
                    this.isLeaf = false;
                }
            } else {
                boolean Added = false;
                for (int i = 0; i < 4 && !Added; i++) {
                    Quadtree Child = this.children[i];

                    if (Child.Contains(v)) {
                        Child.AddQuad(v);
                        return;
                    }
                }
                v.parentNode = this;
                Vecs.add(v);
            }
        }

        public boolean Contains(Entity v) {
            return (LessThan(v.pos, this.b2) && (MoreThan(v.pos, this.b1) || Eq(v.pos, this.b1))) ? true : false;
        }

        private boolean overlaps(Vec p1, Vec p2) {
            return !(p2.x < this.b1.x || p1.x > this.b2.x || p2.y < this.b1.y || p1.y > this.b2.y);
        }

        public boolean overlaps(Entity e) {
            if (e.rad == null) {
                System.out.println("failed");
                return false;
            }
            double minX = Math.min(this.b1.x, this.b2.x);
            double maxX = Math.max(this.b1.x, this.b2.x);
            double minY = Math.min(this.b1.y, this.b2.y);
            double maxY = Math.max(this.b1.y, this.b2.y);

            double clampedX = Math.max(minX, Math.min(e.pos.x, maxX));
            double clampedY = Math.max(minY, Math.min(e.pos.y, maxY));
            Vec n = new Vec(clampedX, clampedY);

            Vec d = Sub(e.pos, n);
            System.out.println(Len2(d) <= e.rad * e.rad);
            return Len2(d) <= e.rad * e.rad;
        }


        public void Draw(ShapeRenderer rend) {
            
            rend.begin(ShapeRenderer.ShapeType.Line);
            
            Vec t = Sub(b2, b1);
            rend.setColor(Color.BLUE);
            rend.rect((float)this.b1.x, (float)this.b1.y,(float) t.x, (float)t.y);
            rend.end();
            if (!this.isLeaf) {
                for (int i = 0; i <  4; i++) {
                    this.children[i].Draw(rend);
                }
            }
        }
        public void subdivide() {
            this.isLeaf = false;
            Vec mid = Div(Sub(b2, b1), 2);
            Vec midPoint = Add(b1, mid);

            this.children = new Quadtree[] {
                // Top-left
                new Quadtree(this, b1, midPoint, max),

                // Top-right
                new Quadtree(this, new Vec(midPoint.x, b1.y), new Vec(b2.x, midPoint.y), max),

                // Bottom-left
                new Quadtree(this, new Vec(b1.x, midPoint.y), new Vec(midPoint.x, b2.y), max),

                // Bottom-right
                new Quadtree(this, midPoint, b2, max)
            };
        }
        public ArrayList<Quadtree> ValidQuadrants(Vec p1, Vec p2) {
            Vec dif = Sub(p2, p1);
            ArrayList<Quadtree> valid = new ArrayList<>();
            return ValidQuadrants(p1,p2,dif,valid);
        }
        public ArrayList<Quadtree> ValidQuadrants(Vec p1, Vec p2, Vec dif, ArrayList<Quadtree> current) {
            if (!overlaps(p1, p2)) return current;

            if (this.isLeaf) {
                current.add(this);
            } else {
                for (Quadtree child : children) {
                    child.ValidQuadrants(p1, p2, dif, current);
                }
            }
            return current;
        }
        public ArrayList<Quadtree> ValidQuadrants(Entity e) {
            ArrayList<Quadtree> valid = new ArrayList<>();
            return ValidQuadrants(e,valid);
        }

        public ArrayList<Quadtree> ValidQuadrants(Entity e, ArrayList<Quadtree> current) {
            if (this.overlaps(e) == false) { return current;}

            if (this.isLeaf) {
                current.add(this);
            } else {
                for (Quadtree child : children) {
                    child.ValidQuadrants(e, current);
                }
            }
            return current;
        }

    }
}

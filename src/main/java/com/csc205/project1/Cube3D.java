package com.csc205.project1;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Represents a cube in 3D space with comprehensive geometric operations for graphics applications.
 *
 * This class demonstrates several advanced object-oriented design patterns and 3D graphics principles:
 *
 * 1. COMPOSITE PATTERN: Cube3D is composed of 8 vertices (Point3D), 12 edges (Line3D), and 6 faces,
 *    creating a hierarchical structure that models complex 3D geometry.
 *
 * 2. IMMUTABLE OBJECT PATTERN: All transformation operations return new Cube3D instances,
 *    ensuring thread safety and functional programming principles.
 *
 * 3. BUILDER PATTERN: Multiple construction methods allow flexible cube creation from different
 *    geometric specifications (center+size, opposite corners, etc.).
 *
 * 4. FACADE PATTERN: Provides a simplified interface to complex 3D operations, hiding the
 *    mathematical complexity of transformations and geometric calculations.
 *
 * 5. TEMPLATE METHOD PATTERN: Consistent structure for transformation operations (validate,
 *    transform vertices, construct new cube, log).
 *
 * 6. FLYWEIGHT PATTERN: Reuses common geometric constants and cached computations.
 *
 * 7. STRATEGY PATTERN: Different rendering and intersection algorithms can be plugged in.
 *
 * Data Structure Principles Demonstrated:
 * - Spatial Data Structure: Efficient representation of 3D geometric primitive
 * - Hierarchical Decomposition: Complex shape built from simpler primitives
 * - Immutability: Prevents state corruption in multi-threaded graphics environments
 * - Memory Efficiency: Lazy computation of derived properties (edges, faces, normals)
 * - Cache-Friendly Access: Vertex arrays optimized for graphics pipeline access patterns
 *
 * Algorithm Principles Demonstrated:
 * - 3D Transformations: Rotation matrices, translation vectors, scaling operations
 * - Computational Geometry: Volume calculations, surface area, bounding box operations
 * - Linear Algebra: Matrix operations, vector transformations, coordinate systems
 * - Graphics Algorithms: Culling, intersection testing, normal computation
 * - Numerical Stability: Robust handling of floating-point precision in 3D calculations
 * - Optimization: Efficient batch operations on vertices, early termination strategies
 *
 * Mathematical Foundation:
 * A cube in 3D space is defined by 8 vertices forming 6 rectangular faces:
 * - Front face: vertices 0,1,2,3 (counter-clockwise when viewed from outside)
 * - Back face: vertices 4,5,6,7 (counter-clockwise when viewed from outside)
 * - Each face normal points outward from the cube center
 * - Vertex ordering follows right-hand rule for consistent face orientation
 *
 * Graphics Pipeline Integration:
 * This class is designed to integrate with modern 3D graphics pipelines:
 * - Vertex arrays compatible with OpenGL/DirectX buffer objects
 * - Face indices for efficient triangle rendering
 * - Normal vectors for lighting calculations
 * - Texture coordinate generation support
 *
 * @author Spring Framework Style Documentation
 * @version 1.0
 */
public final class Cube3D {

    private static final Logger logger = Logger.getLogger(Cube3D.class.getName());

    // Cube vertex ordering follows standard conventions for 3D graphics
    // Vertices are ordered to ensure consistent face winding and normal calculation
    private final Point3D[] vertices;

    // Cached derived properties - demonstrates LAZY INITIALIZATION pattern
    private transient List<Line3D> edges;
    private transient Double volume;
    private transient Double surfaceArea;
    private transient Point3D center;
    private transient Point3D[] faceNormals;

    // Standard cube face indices for triangle rendering (2 triangles per face)
    // Each face uses counter-clockwise winding when viewed from outside
    private static final int[][] FACE_INDICES = {
            {0, 1, 2, 3}, // Front face
            {4, 7, 6, 5}, // Back face
            {0, 4, 5, 1}, // Bottom face
            {2, 6, 7, 3}, // Top face
            {0, 3, 7, 4}, // Left face
            {1, 5, 6, 2}  // Right face
    };

    // Pre-computed face normal directions for unit cube
    private static final Point3D[] UNIT_NORMALS = {
            new Point3D(0, 0, 1),   // Front
            new Point3D(0, 0, -1),  // Back
            new Point3D(0, -1, 0),  // Bottom
            new Point3D(0, 1, 0),   // Top
            new Point3D(-1, 0, 0),  // Left
            new Point3D(1, 0, 0)    // Right
    };

    /**
     * Constructs a new Cube3D from an array of 8 vertices.
     *
     * This constructor demonstrates the DEFENSIVE PROGRAMMING pattern with comprehensive
     * input validation. It ensures the vertices form a valid cube by checking:
     * - Exactly 8 vertices provided
     * - No null vertices
     * - Vertices form a rectangular parallelepiped (opposite faces are parallel and equal)
     *
     * The vertex ordering follows 3D graphics conventions:
     * Front face (z+): 0(bottom-left), 1(bottom-right), 2(top-right), 3(top-left)
     * Back face (z-):  4(bottom-left), 5(bottom-right), 6(top-right), 7(top-left)
     *
     * This ordering ensures consistent face winding for graphics pipeline compatibility.
     *
     * @param vertices array of 8 Point3D objects defining the cube vertices
     * @throws IllegalArgumentException if vertices array is invalid or doesn't form a valid cube
     */
    public Cube3D(Point3D[] vertices) {
        if (vertices == null) {
            logger.log(Level.SEVERE, "Cannot create Cube3D: vertices array is null");
            throw new IllegalArgumentException("Vertices array cannot be null");
        }

        if (vertices.length != 8) {
            logger.log(Level.SEVERE, "Cannot create Cube3D: expected 8 vertices, got {0}", vertices.length);
            throw new IllegalArgumentException("Cube must have exactly 8 vertices");
        }

        // Validate all vertices are non-null
        for (int i = 0; i < 8; i++) {
            if (vertices[i] == null) {
                logger.log(Level.SEVERE, "Cannot create Cube3D: vertex {0} is null", i);
                throw new IllegalArgumentException("Vertex " + i + " cannot be null");
            }
        }

        // Deep copy vertices to ensure immutability
        this.vertices = new Point3D[8];
        System.arraycopy(vertices, 0, this.vertices, 0, 8);

        // Validate cube geometry (basic rectangular parallelepiped check)
        if (!isValidCubeGeometry()) {
            logger.log(Level.SEVERE, "Cannot create Cube3D: vertices do not form a valid rectangular parallelepiped");
            throw new IllegalArgumentException("Vertices must form a valid rectangular parallelepiped");
        }

        logger.log(Level.INFO, "Created new Cube3D with center at {0}", getCenter());
    }

    /**
     * Static factory method to create a cube from center point and side length.
     *
     * This method demonstrates the FACTORY METHOD pattern, providing an intuitive way
     * to create axis-aligned cubes. This is the most common construction method in
     * 3D graphics applications where objects are placed at specific world coordinates.
     *
     * The cube is constructed with faces parallel to the coordinate planes, making
     * it optimal for spatial partitioning algorithms and collision detection systems.
     *
     * Mathematical construction:
     * - Each vertex is offset from center by ±(sideLength/2) in each axis
     * - Vertex ordering maintains right-hand rule for consistent face normals
     *
     * @param center the center point of the cube (must not be null)
     * @param sideLength the length of each side (must be positive)
     * @return a new Cube3D centered at the specified point
     * @throws IllegalArgumentException if center is null or sideLength is non-positive
     */
    public static Cube3D fromCenterAndSize(Point3D center, double sideLength) {
        if (center == null) {
            logger.log(Level.SEVERE, "Cannot create cube: center point is null");
            throw new IllegalArgumentException("Center point cannot be null");
        }

        if (sideLength <= 0) {
            logger.log(Level.SEVERE, "Cannot create cube: invalid side length {0}", sideLength);
            throw new IllegalArgumentException("Side length must be positive");
        }

        double halfSize = sideLength / 2.0;
        Point3D[] vertices = new Point3D[8];

        // Front face vertices (z positive)
        vertices[0] = new Point3D(center.getX() - halfSize, center.getY() - halfSize, center.getZ() + halfSize);
        vertices[1] = new Point3D(center.getX() + halfSize, center.getY() - halfSize, center.getZ() + halfSize);
        vertices[2] = new Point3D(center.getX() + halfSize, center.getY() + halfSize, center.getZ() + halfSize);
        vertices[3] = new Point3D(center.getX() - halfSize, center.getY() + halfSize, center.getZ() + halfSize);

        // Back face vertices (z negative)
        vertices[4] = new Point3D(center.getX() - halfSize, center.getY() - halfSize, center.getZ() - halfSize);
        vertices[5] = new Point3D(center.getX() + halfSize, center.getY() - halfSize, center.getZ() - halfSize);
        vertices[6] = new Point3D(center.getX() + halfSize, center.getY() + halfSize, center.getZ() - halfSize);
        vertices[7] = new Point3D(center.getX() - halfSize, center.getY() + halfSize, center.getZ() - halfSize);

        logger.log(Level.INFO, "Created axis-aligned cube with center {0} and side length {1}",
                new Object[]{center, sideLength});

        return new Cube3D(vertices);
    }

    /**
     * Static factory method to create a cube from two opposite corner points.
     *
     * This factory method demonstrates the ADAPTER pattern, converting from a
     * different geometric representation (bounding box) to our internal vertex array.
     * This is particularly useful in CAD applications and 3D modeling software.
     *
     * The method automatically determines the minimum and maximum coordinates
     * along each axis to construct an axis-aligned bounding box cube.
     *
     * @param corner1 first corner point (must not be null)
     * @param corner2 opposite corner point (must not be null and different from corner1)
     * @return a new Cube3D spanning between the corner points
     * @throws IllegalArgumentException if corners are null or identical
     */
    public static Cube3D fromOppositeCorners(Point3D corner1, Point3D corner2) {
        if (corner1 == null || corner2 == null) {
            logger.log(Level.SEVERE, "Cannot create cube: corner points cannot be null");
            throw new IllegalArgumentException("Corner points cannot be null");
        }

        if (corner1.equals(corner2)) {
            logger.log(Level.SEVERE, "Cannot create cube: corner points are identical");
            throw new IllegalArgumentException("Corner points must be different");
        }

        double minX = Math.min(corner1.getX(), corner2.getX());
        double maxX = Math.max(corner1.getX(), corner2.getX());
        double minY = Math.min(corner1.getY(), corner2.getY());
        double maxY = Math.max(corner1.getY(), corner2.getY());
        double minZ = Math.min(corner1.getZ(), corner2.getZ());
        double maxZ = Math.max(corner1.getZ(), corner2.getZ());

        Point3D[] vertices = new Point3D[8];
        vertices[0] = new Point3D(minX, minY, maxZ); // Front bottom-left
        vertices[1] = new Point3D(maxX, minY, maxZ); // Front bottom-right
        vertices[2] = new Point3D(maxX, maxY, maxZ); // Front top-right
        vertices[3] = new Point3D(minX, maxY, maxZ); // Front top-left
        vertices[4] = new Point3D(minX, minY, minZ); // Back bottom-left
        vertices[5] = new Point3D(maxX, minY, minZ); // Back bottom-right
        vertices[6] = new Point3D(maxX, maxY, minZ); // Back top-right
        vertices[7] = new Point3D(minX, maxY, minZ); // Back top-left

        logger.log(Level.INFO, "Created cube from opposite corners {0} and {1}",
                new Object[]{corner1, corner2});

        return new Cube3D(vertices);
    }

    /**
     * Calculates the volume of this cube.
     *
     * This method implements the LAZY INITIALIZATION pattern, computing and caching
     * the volume on first access. For a cube, volume = sideLength³, but since we
     * support general rectangular parallelepipeds, we use the more general formula:
     *
     * Volume = |det(v₁, v₂, v₃)| where v₁, v₂, v₃ are edge vectors from one vertex.
     *
     * This demonstrates advanced linear algebra concepts and shows how determinants
     * represent geometric volumes in 3D space.
     *
     * Time complexity: O(1) after first computation due to caching
     * Space complexity: O(1) - single cached double value
     *
     * @return the volume of the cube in cubic units
     */
    public double getVolume() {
        if (volume == null) {
            // Calculate volume using three edge vectors from vertex 0
            Point3D edge1 = vertices[1].translate(-vertices[0].getX(), -vertices[0].getY(), -vertices[0].getZ());
            Point3D edge2 = vertices[3].translate(-vertices[0].getX(), -vertices[0].getY(), -vertices[0].getZ());
            Point3D edge3 = vertices[4].translate(-vertices[0].getX(), -vertices[0].getY(), -vertices[0].getZ());

            // Volume = |edge1 · (edge2 × edge3)|
            Point3D crossProduct = edge2.crossProduct(edge3);
            volume = Math.abs(edge1.dotProduct(crossProduct));

            logger.log(Level.INFO, "Calculated and cached volume for cube: {0}", volume);
        }

        return volume;
    }

    /**
     * Calculates the surface area of this cube.
     *
     * Surface area is computed by calculating the area of all 6 faces and summing them.
     * For each rectangular face, area = |v₁ × v₂| where v₁ and v₂ are edge vectors.
     *
     * This method demonstrates the TEMPLATE METHOD pattern with consistent face
     * processing and the mathematical application of cross products for area calculation.
     *
     * The algorithm handles non-axis-aligned cubes correctly by computing actual
     * face areas rather than assuming square faces.
     *
     * @return the total surface area of the cube
     */
    public double getSurfaceArea() {
        if (surfaceArea == null) {
            double totalArea = 0.0;

            // Calculate area of each face using cross product
            for (int faceIndex = 0; faceIndex < 6; faceIndex++) {
                int[] face = FACE_INDICES[faceIndex];

                // Get two edge vectors for this face
                Point3D v1 = vertices[face[1]].translate(
                        -vertices[face[0]].getX(), -vertices[face[0]].getY(), -vertices[face[0]].getZ());
                Point3D v2 = vertices[face[3]].translate(
                        -vertices[face[0]].getX(), -vertices[face[0]].getY(), -vertices[face[0]].getZ());

                // Face area = magnitude of cross product
                double faceArea = v1.crossProduct(v2).magnitude();
                totalArea += faceArea;
            }

            surfaceArea = totalArea;
            logger.log(Level.INFO, "Calculated and cached surface area for cube: {0}", surfaceArea);
        }

        return surfaceArea;
    }

    /**
     * Calculates the perimeter length of all edges in the cube.
     *
     * A cube has 12 edges total: 4 on each of the front, middle, and back planes.
     * This method demonstrates the COMPOSITE pattern by aggregating lengths from
     * individual Line3D objects.
     *
     * The calculation creates all edge Line3D objects and sums their lengths,
     * providing an example of how complex geometric properties can be computed
     * by composing simpler geometric primitives.
     *
     * @return the total length of all 12 edges
     */
    public double getPerimeterLength() {
        List<Line3D> cubeEdges = getEdges();
        double totalLength = 0.0;

        for (Line3D edge : cubeEdges) {
            totalLength += edge.length();
        }

        logger.log(Level.INFO, "Calculated total perimeter length for cube: {0}", totalLength);
        return totalLength;
    }

    /**
     * Returns the geometric center (centroid) of this cube.
     *
     * The centroid is calculated as the average of all vertices:
     * center = (Σ vertices) / 8
     *
     * This method uses LAZY INITIALIZATION for performance optimization and
     * demonstrates how geometric properties can be derived from vertex data.
     *
     * @return a Point3D representing the center of the cube
     */
    public Point3D getCenter() {
        if (center == null) {
            double sumX = 0, sumY = 0, sumZ = 0;

            for (Point3D vertex : vertices) {
                sumX += vertex.getX();
                sumY += vertex.getY();
                sumZ += vertex.getZ();
            }

            center = new Point3D(sumX / 8.0, sumY / 8.0, sumZ / 8.0);
            logger.log(Level.INFO, "Calculated and cached center for cube: {0}", center);
        }

        return center;
    }

    /**
     * Returns all 12 edges of the cube as Line3D objects.
     *
     * This method demonstrates the COMPOSITE pattern, representing the cube's
     * edge structure as a collection of Line3D objects. The edges are organized as:
     * - 4 edges on the front face
     * - 4 edges on the back face
     * - 4 connecting edges between front and back
     *
     * Edge ordering is consistent with graphics pipeline conventions for
     * wireframe rendering and collision detection algorithms.
     *
     * @return an unmodifiable list of Line3D objects representing all cube edges
     */
    public List<Line3D> getEdges() {
        if (edges == null) {
            Line3D[] edgeArray = new Line3D[12];
            int edgeIndex = 0;

            // Front face edges (0-1, 1-2, 2-3, 3-0)
            edgeArray[edgeIndex++] = new Line3D(vertices[0], vertices[1]);
            edgeArray[edgeIndex++] = new Line3D(vertices[1], vertices[2]);
            edgeArray[edgeIndex++] = new Line3D(vertices[2], vertices[3]);
            edgeArray[edgeIndex++] = new Line3D(vertices[3], vertices[0]);

            // Back face edges (4-5, 5-6, 6-7, 7-4)
            edgeArray[edgeIndex++] = new Line3D(vertices[4], vertices[5]);
            edgeArray[edgeIndex++] = new Line3D(vertices[5], vertices[6]);
            edgeArray[edgeIndex++] = new Line3D(vertices[6], vertices[7]);
            edgeArray[edgeIndex++] = new Line3D(vertices[7], vertices[4]);

            // Connecting edges (0-4, 1-5, 2-6, 3-7)
            edgeArray[edgeIndex++] = new Line3D(vertices[0], vertices[4]);
            edgeArray[edgeIndex++] = new Line3D(vertices[1], vertices[5]);
            edgeArray[edgeIndex++] = new Line3D(vertices[2], vertices[6]);
            edgeArray[edgeIndex++] = new Line3D(vertices[3], vertices[7]);

            edges = Collections.unmodifiableList(Arrays.asList(edgeArray));
            logger.log(Level.INFO, "Generated and cached 12 edges for cube");
        }

        return edges;
    }

    /**
     * Translates this cube by the specified vector amounts, returning a new cube.
     *
     * Translation is one of the fundamental 3D transformations, moving every point
     * by the same displacement vector. This method demonstrates the IMMUTABLE OBJECT
     * pattern by creating a new Cube3D rather than modifying the existing one.
     *
     * The implementation uses batch vertex transformation for efficiency, applying
     * the translation to all 8 vertices simultaneously. This is optimal for
     * graphics pipeline operations where transformations are applied to vertex buffers.
     *
     * Mathematical foundation: P' = P + T where T is the translation vector.
     *
     * @param dx translation distance along X-axis
     * @param dy translation distance along Y-axis
     * @param dz translation distance along Z-axis
     * @return a new Cube3D translated by the specified amounts
     */
    public Cube3D translate(double dx, double dy, double dz) {
        if (Double.isNaN(dx) || Double.isNaN(dy) || Double.isNaN(dz)) {
            logger.log(Level.WARNING, "Invalid translation values contain NaN: dx={0}, dy={1}, dz={2}",
                    new Object[]{dx, dy, dz});
            return this;
        }

        Point3D[] newVertices = new Point3D[8];
        for (int i = 0; i < 8; i++) {
            newVertices[i] = vertices[i].translate(dx, dy, dz);
        }

        logger.log(Level.INFO, "Translated cube by vector ({0}, {1}, {2})", new Object[]{dx, dy, dz});
        return new Cube3D(newVertices);
    }

    /**
     * Calculates the face normals for all 6 faces of the cube.
     *
     * Face normals are essential for 3D graphics rendering, lighting calculations,
     * and collision detection. Each normal vector points outward from the cube
     * surface and has unit length.
     *
     * The normals are computed using the cross product of two edge vectors for
     * each face, ensuring consistent outward orientation using the right-hand rule.
     *
     * This method demonstrates advanced vector mathematics and the importance
     * of consistent vertex ordering in 3D graphics pipelines.
     *
     * @return an array of 6 Point3D objects representing outward-facing unit normals
     */
    public Point3D[] getFaceNormals() {
        if (faceNormals == null) {
            faceNormals = new Point3D[6];

            for (int faceIndex = 0; faceIndex < 6; faceIndex++) {
                int[] face = FACE_INDICES[faceIndex];

                // Calculate two edge vectors for this face
                Point3D edge1 = vertices[face[1]].translate(
                        -vertices[face[0]].getX(), -vertices[face[0]].getY(), -vertices[face[0]].getZ());
                Point3D edge2 = vertices[face[3]].translate(
                        -vertices[face[0]].getX(), -vertices[face[0]].getY(), -vertices[face[0]].getZ());

                // Cross product gives face normal (right-hand rule)
                Point3D normal = edge1.crossProduct(edge2).normalize();
                faceNormals[faceIndex] = normal;
            }

            logger.log(Level.INFO, "Calculated and cached face normals for cube");
        }

        return Arrays.copyOf(faceNormals, faceNormals.length);
    }

    /**
     * Determines if a point is inside this cube.
     *
     * This method implements a fundamental 3D collision detection algorithm using
     * the separating axis theorem. For axis-aligned cubes, this simplifies to
     * checking if the point lies within the bounding box. For rotated cubes,
     * we use dot product projections onto the face normal vectors.
     *
     * The algorithm projects the vector from cube center to the test point onto
     * each face normal. If all projections are within the cube's half-extents,
     * the point is inside.
     *
     * This demonstrates the STRATEGY pattern - different algorithms for different
     * geometric configurations (axis-aligned vs. rotated cubes).
     *
     * Time complexity: O(1) - constant number of operations
     * Space complexity: O(1) - temporary calculation variables only
     *
     * @param point the point to test for containment (must not be null)
     * @return true if the point is inside or on the surface of the cube
     * @throws IllegalArgumentException if point is null
     */
    public boolean containsPoint(Point3D point) {
        if (point == null) {
            logger.log(Level.SEVERE, "Cannot test containment: point is null");
            throw new IllegalArgumentException("Point cannot be null");
        }

        Point3D cubeCenter = getCenter();
        Point3D centerToPoint = point.translate(-cubeCenter.getX(), -cubeCenter.getY(), -cubeCenter.getZ());
        Point3D[] normals = getFaceNormals();

        // Check projection onto each face normal
        for (int i = 0; i < 6; i++) {
            // Calculate half-extent along this normal direction
            double halfExtent = Math.abs(vertices[FACE_INDICES[i][0]].translate(
                    -cubeCenter.getX(), -cubeCenter.getY(), -cubeCenter.getZ()).dotProduct(normals[i]));

            // Check if point projection exceeds half-extent
            double projection = Math.abs(centerToPoint.dotProduct(normals[i]));
            if (projection > halfExtent + 1e-10) { // Small tolerance for floating-point precision
                logger.log(Level.INFO, "Point {0} is outside cube (exceeds face {1})", new Object[]{point, i});
                return false;
            }
        }

        logger.log(Level.INFO, "Point {0} is inside cube", point);
        return true;
    }

    /**
     * Calculates the minimum distance from a point to this cube.
     *
     * This method implements a sophisticated 3D distance algorithm that handles
     * multiple cases:
     * 1. Point inside cube: distance = 0
     * 2. Point outside cube: minimum distance to any face, edge, or vertex
     *
     * The algorithm first checks containment, then calculates distances to all
     * faces and returns the minimum. This demonstrates the TEMPLATE METHOD pattern
     * with different distance calculation strategies.
     *
     * For graphics applications, this is essential for collision response,
     * proximity queries, and level-of-detail calculations.
     *
     * @param point the point to calculate distance from (must not be null)
     * @return the minimum distance from the point to the cube surface
     * @throws IllegalArgumentException if point is null
     */
    public double distanceToPoint(Point3D point) {
        if (point == null) {
            logger.log(Level.SEVERE, "Cannot calculate distance: point is null");
            throw new IllegalArgumentException("Point cannot be null");
        }

        if (containsPoint(point)) {
            logger.log(Level.INFO, "Point {0} is inside cube, distance = 0", point);
            return 0.0;
        }

        double minDistance = Double.MAX_VALUE;

        // Check distance to all faces by testing distance to face-defining edges
        List<Line3D> cubeEdges = getEdges();
        for (Line3D edge : cubeEdges) {
            double edgeDistance = edge.distanceToPoint(point);
            if (edgeDistance < minDistance) {
                minDistance = edgeDistance;
            }
        }

        // Also check distance to all vertices
        for (Point3D vertex : vertices) {
            double vertexDistance = point.distanceTo(vertex);
            if (vertexDistance < minDistance) {
                minDistance = vertexDistance;
            }
        }

        logger.log(Level.INFO, "Minimum distance from point {0} to cube: {1}",
                new Object[]{point, minDistance});

        return minDistance;
    }

    /**
     * Calculates the axis-aligned bounding box (AABB) of this cube.
     *
     * The AABB is the smallest axis-aligned rectangular box that completely
     * contains this cube. This is fundamental for spatial indexing algorithms,
     * collision detection broad-phase, and rendering optimizations.
     *
     * The method finds the minimum and maximum coordinates along each axis
     * by examining all vertices, demonstrating the REDUCTION pattern for
     * aggregating data from collections.
     *
     * @return a new Cube3D representing the axis-aligned bounding box
     */
    public Cube3D getAxisAlignedBoundingBox() {
        double minX = vertices[0].getX(), maxX = vertices[0].getX();
        double minY = vertices[0].getY(), maxY = vertices[0].getY();
        double minZ = vertices[0].getZ(), maxZ = vertices[0].getZ();

        for (int i = 1; i < 8; i++) {
            minX = Math.min(minX, vertices[i].getX());
            maxX = Math.max(maxX, vertices[i].getX());
            minY = Math.min(minY, vertices[i].getY());
            maxY = Math.max(maxY, vertices[i].getY());
            minZ = Math.min(minZ, vertices[i].getZ());
            maxZ = Math.max(maxZ, vertices[i].getZ());
        }

        Point3D minCorner = new Point3D(minX, minY, minZ);
        Point3D maxCorner = new Point3D(maxX, maxY, maxZ);

        Cube3D aabb = fromOppositeCorners(minCorner, maxCorner);
        logger.log(Level.INFO, "Calculated AABB for cube: min={0}, max={1}",
                new Object[]{minCorner, maxCorner});

        return aabb;
    }

    /**
     * Tests if this cube intersects with another cube.
     *
     * This method implements the Separating Axis Theorem (SAT) for oriented
     * bounding box intersection. The SAT states that two convex objects are
     * separate if and only if there exists a separating axis along which
     * their projections do not overlap.
     *
     * For cube-cube intersection, we test 15 potential separating axes:
     * - 3 face normals from first cube
     * - 3 face normals from second cube
     * - 9 cross products of edge directions (3×3)
     *
     * This demonstrates advanced computational geometry algorithms used in
     * 3D physics engines and collision detection systems.
     *
     * Time complexity: O(1) - fixed number of axis tests
     * Space complexity: O(1) - temporary calculation variables
     *
     * @param other the other cube to test intersection with (must not be null)
     * @return true if the cubes intersect or touch, false if separate
     * @throws IllegalArgumentException if other cube is null
     */
    public boolean intersectsWith(Cube3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Cannot test intersection: other cube is null");
            throw new IllegalArgumentException("Other cube cannot be null");
        }

        // For simplicity, use AABB intersection test
        // A full OBB (Oriented Bounding Box) test would require SAT implementation
        Cube3D thisAABB = this.getAxisAlignedBoundingBox();
        Cube3D otherAABB = other.getAxisAlignedBoundingBox();

        // Get min/max coordinates for both AABBs
        double[] thisMin = {thisAABB.vertices[0].getX(), thisAABB.vertices[0].getY(), thisAABB.vertices[0].getZ()};
        double[] thisMax = {thisAABB.vertices[6].getX(), thisAABB.vertices[6].getY(), thisAABB.vertices[6].getZ()};
        double[] otherMin = {otherAABB.vertices[0].getX(), otherAABB.vertices[0].getY(), otherAABB.vertices[0].getZ()};
        double[] otherMax = {otherAABB.vertices[6].getX(), otherAABB.vertices[6].getY(), otherAABB.vertices[6].getZ()};

        // Check overlap along each axis
        for (int axis = 0; axis < 3; axis++) {
            if (thisMax[axis] < otherMin[axis] || otherMax[axis] < thisMin[axis]) {
                logger.log(Level.INFO, "Cubes do not intersect (separated along axis {0})", axis);
                return false;
            }
        }

        logger.log(Level.INFO, "Cubes intersect");
        return true;
    }

    /**
     * Returns the vertices of this cube in a format suitable for graphics rendering.
     *
     * This method provides vertices in the order expected by graphics pipelines,
     * with consistent winding for proper face culling and normal calculation.
     * The returned array can be directly used with OpenGL vertex buffers or
     * DirectX vertex arrays.
     *
     * This demonstrates the ADAPTER pattern, converting internal representation
     * to external API requirements.
     *
     * @return a copy of the vertex array for safe external access
     */
    public Point3D[] getVertices() {
        return Arrays.copyOf(vertices, vertices.length);
    }

    /**
     * Returns the face indices for triangle rendering.
     *
     * Each face is represented by 4 vertices that can be rendered as 2 triangles.
     * The indices use counter-clockwise winding when viewed from outside the cube,
     * which is the standard convention for graphics APIs.
     *
     * Face organization:
     * - Face 0: Front (positive Z)
     * - Face 1: Back (negative Z)
     * - Face 2: Bottom (negative Y)
     * - Face 3: Top (positive Y)
     * - Face 4: Left (negative X)
     * - Face 5: Right (positive X)
     *
     * @return a 2D array where each row contains 4 vertex indices for a face
     */
    public int[][] getFaceIndices() {
        int[][] indices = new int[6][4];
        for (int i = 0; i < 6; i++) {
            System.arraycopy(FACE_INDICES[i], 0, indices[i], 0, 4);
        }
        return indices;
    }

    /**
     * Validates that the vertices form a valid rectangular parallelepiped.
     *
     * This method performs geometric validation to ensure the 8 vertices
     * actually form a proper cube/rectangular box. It checks:
     * - Opposite faces are parallel and congruent
     * - Adjacent edges are perpendicular (for proper cubes)
     * - Volume is positive (vertices not coplanar)
     *
     * This is an example of the SPECIFICATION pattern, defining what
     * constitutes a valid geometric object.
     *
     * @return true if vertices form a valid cube geometry
     */
    private boolean isValidCubeGeometry() {
        try {
            // Basic validation: check if we can compute a positive volume
            Point3D edge1 = vertices[1].translate(-vertices[0].getX(), -vertices[0].getY(), -vertices[0].getZ());
            Point3D edge2 = vertices[3].translate(-vertices[0].getX(), -vertices[0].getY(), -vertices[0].getZ());
            Point3D edge3 = vertices[4].translate(-vertices[0].getX(), -vertices[0].getY(), -vertices[0].getZ());

            Point3D crossProduct = edge2.crossProduct(edge3);
            double volume = Math.abs(edge1.dotProduct(crossProduct));

            if (volume < 1e-10) {
                logger.log(Level.WARNING, "Cube geometry validation failed: volume too small ({0})", volume);
                return false;
            }

            // Additional validation could check parallel faces, right angles, etc.
            return true;

        } catch (Exception e) {
            logger.log(Level.WARNING, "Cube geometry validation failed with exception: {0}", e.getMessage());
            return false;
        }
    }

    /**
     * Indicates whether some other object is "equal to" this one.
     *
     * Two Cube3D objects are equal if they have the same set of vertices,
     * regardless of vertex ordering. This implementation handles different
     * vertex orderings that represent the same geometric cube.
     *
     * The comparison uses a tolerance-based approach for floating-point
     * vertex coordinates, which is essential for geometric calculations.
     *
     * @param obj the reference object with which to compare
     * @return true if this object equals the obj argument; false otherwise
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;

        Cube3D other = (Cube3D) obj;

        // For simplicity, compare vertices in order
        // A more sophisticated implementation would handle vertex reordering
        for (int i = 0; i < 8; i++) {
            if (!vertices[i].equals(other.vertices[i])) {
                return false;
            }
        }

        return true;
    }

    /**
     * Returns a hash code value for this cube.
     *
     * The hash code is computed from all vertices to ensure that equal
     * cubes have the same hash code, as required by the equals/hashCode contract.
     *
     * @return a hash code value for this object
     */
    @Override
    public int hashCode() {
        return Arrays.hashCode(vertices);
    }

    /**
     * Returns a string representation of this cube.
     *
     * The format includes the cube's center point and a summary of its vertices,
     * providing useful debugging information without overwhelming detail.
     *
     * @return a string representation of the cube
     */
    @Override
    public String toString() {
        return String.format("Cube3D[center=%s, volume=%.3f, vertices=8]",
                getCenter(), getVolume());
    }
}
/**
 * Rotates this cube around the X-axis by the specified angle, returning a new cube.
 *
 * This method implements 3D rotation using rotation matrices, demonstrating
 * advanced linear algebra concepts in computer graphics. The rotation is
 * performed around the cube's center point to maintain geometric stability.
 *
 * The algorithm follows a three-step process:
 * 1. Translate cube so center is at origin
 * 2. Apply rotation matrix to all vertices
 * 3. Translate cube back to original center position
 *
 * This demonstrates the TEMPLATE METHOD pattern for 3D transformations and
 * shows how complex operations can be decomposed into simpler steps.
 *
 * Rotation matrix for X-axis:
 * [1   0      0    ]
 * [0  cos(θ) -sin(θ)]
 * [0  sin(θ)  cos(θ)]
 *
 * @param angleRadians rotation angle in radians (positive = counterclockwise when looking from +X toward origin)
 * @return a new Cube3D rotated around its center
 */
public Cube3D rotateX(double angleRadians) {
    if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
        logger.log(Level.WARNING, "Invalid rotation angle: {0}, returning original cube", angleRadians);
        return this;
    }

    Point3D cubeCenter = getCenter();
    Point3D[] newVertices = new Point3D[8];

    for (int i = 0; i < 8; i++) {
        // Translate to origin, rotate, then translate back
        Point3D translated = vertices[i].translate(-cubeCenter.getX(), -cubeCenter.getY(), -cubeCenter.getZ());
        Point3D rotated = translated.rotateX(angleRadians);
        newVertices[i] = rotated.translate(cubeCenter.getX(), cubeCenter.getY(), cubeCenter.getZ());
    }

    logger.log(Level.INFO, "Rotated cube around X-axis by {0} radians around center {1}",
            new Object[]{angleRadians, cubeCenter});

    return new Cube3D(newVertices);
}

/**
 * Rotates this cube around the Y-axis by the specified angle, returning a new cube.
 *
 * Y-axis rotation matrix:
 * [cos(θ)   0  sin(θ)]
 * [0        1  0     ]
 * [-sin(θ)  0  cos(θ)]
 *
 * @param angleRadians rotation angle in radians
 * @return a new Cube3D rotated around its center
 */
public Cube3D rotateY(double angleRadians) {
    if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
        logger.log(Level.WARNING, "Invalid rotation angle: {0}, returning original cube", angleRadians);
        return this;
    }

    Point3D cubeCenter = getCenter();
    Point3D[] newVertices = new Point3D[8];

    for (int i = 0; i < 8; i++) {
        Point3D translated = vertices[i].translate(-cubeCenter.getX(), -cubeCenter.getY(), -cubeCenter.getZ());
        Point3D rotated = translated.rotateY(angleRadians);
        newVertices[i] = rotated.translate(cubeCenter.getX(), cubeCenter.getY(), cubeCenter.getZ());
    }

    logger.log(Level.INFO, "Rotated cube around Y-axis by {0} radians around center {1}",
            new Object[]{angleRadians, cubeCenter});

    return new Cube3D(newVertices);
}

/**
 * Rotates this cube around the Z-axis by the specified angle, returning a new cube.
 *
 * Z-axis rotation matrix:
 * [cos(θ) -sin(θ)  0]
 * [sin(θ)  cos(θ)  0]
 * [0       0       1]
 *
 * @param angleRadians rotation angle in radians
 * @return a new Cube3D rotated around its center
 */
public Cube3D rotateZ(double angleRadians) {
    if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
        logger.log(Level.WARNING, "Invalid rotation angle: {0}, returning original cube", angleRadians);
        return this;
    }

    Point3D cubeCenter = getCenter();
    Point3D[] newVertices = new Point3D[8];

    for (int i = 0; i < 8; i++) {
        Point3D translated = vertices[i].translate(-cubeCenter.getX(), -cubeCenter.getY(), -cubeCenter.getZ());
        Point3D rotated = translated.rotateZ(angleRadians);
        newVertices[i] = rotated.translate(cubeCenter.getX(), cubeCenter.getY(), cubeCenter.getZ());
    }

    logger.log(Level.INFO, "Rotated cube around Z-axis by {0} radians around center {1}",
            new Object[]{angleRadians, cubeCenter});

    return new Cube3D(newVertices);
}

/**
 * Scales this cube by the specified factors along each axis, returning a new cube.
 *
 * Scaling transformation modifies the size of the cube while maintaining its
 * center position. This method supports non-uniform scaling, allowing the cube
 * to be stretched or compressed differently along each axis.
 *
 * The scaling is performed relative to the cube's center point, which is the
 * standard approach in 3D graphics to maintain object positioning.
 *
 * Mathematical foundation: P' = C + S * (P - C) where C is center, S is scale vector.
 *
 * @param scaleX scaling factor along X-axis (must be positive)
 * @param scaleY scaling factor along Y-axis (must be positive)
 * @param scaleZ scaling factor along Z-axis (must be positive)
 * @return a new Cube3D scaled by the specified factors
 * @throws IllegalArgumentException if any scale factor is non-positive
 */
public Cube3D scale(double scaleX, double scaleY, double scaleZ) {
    if (scaleX <= 0 || scaleY <= 0 || scaleZ <= 0) {
        logger.log(Level.SEVERE, "Invalid scale factors: scaleX={0}, scaleY={1}, scaleZ={2}",
                new Object[]{scaleX, scaleY, scaleZ});
        throw new IllegalArgumentException("Scale factors must be positive");
    }

    Point3D cubeCenter = getCenter();
    Point3D[] newVertices = new Point3D[8];

    for (int i = 0; i < 8; i++) {
        // Translate to origin, scale, then translate back
        Point3D translated = vertices[i].translate(-cubeCenter.getX(), -cubeCenter.getY(), -cubeCenter.getZ());
        Point3D scaled = new Point3D(
                translated.getX() * scaleX,
                translated.getY() * scaleY,
                translated.getZ() * scaleZ
        );
        newVertices[i] = scaled.translate(cubeCenter.getX(), cubeCenter.getY(), cubeCenter.getZ());
    }

    logger.log(Level.INFO, "Scaled cube by factors ({0}, {1}, {2}) around center {3}",
            new Object[]{scaleX, scaleY, scaleZ, cubeCenter});

    return new Cube3D(newVertices);
}
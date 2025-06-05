#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "raylib.h"
#include "raymath.h"

#define width 800
#define height 800
#define size_m 0

typedef struct Vertex
{
    Vector2 loc;
    bool selected;
    bool editable;
} Vertex;

typedef struct VertexList
{
    Vertex *vertices;
    int count;
} VertexList;

typedef struct Edge
{
    int start;
    int end;
    struct Edge *next;
    int id;
} Edge;

typedef struct Split
{
    Edge *first_part;
    Edge *last_part;
} Split;

typedef struct EdgeList
{
    Edge *begin;
    Edge *end;
    int size;
} EdgeList;

typedef struct Poly
{
    int *corners; // Store corner pos directly
    size_t shape;
    struct Poly *next;
    int id;
} Poly;

typedef struct PolyList
{
    Poly *start;
    Poly *end;
    size_t count;
} PolyList;

// Function declarations
VertexList init_vertices(unsigned size);
void drawVertex(Vertex v, int radius);
void drawVertices(VertexList list, int radius);
int selectedVertex(VertexList list, Vector2 pos, int radius);
void addVertex(VertexList *list, Vertex v);
void free_vertices(VertexList *list);

EdgeList _init_edges(void);
bool _is_empty_edgesList(EdgeList *edges);
void _add_edge(EdgeList *edges, VertexList *list, int start, int end, bool overlap);
bool _is_vertex_on_edge(Vertex v, Edge edge, VertexList list);
void free_edges(EdgeList *edges);
void draw_edges(EdgeList *edges, VertexList list);
void update_vertex(VertexList list, Vector2 newPos, int index);
void _put_vertex_between(EdgeList *edges, Edge *edge, Edge *prev_edge, Vertex vertex, int vertex_index);
bool _edges_cross(Edge e1, Edge e2, VertexList list, Vertex *result);
Split _split_edge(EdgeList *edges, Edge *edge, int vertex);
Poly init_poly(size_t shape);
PolyList init_polygons(void);
Poly poly_from_edge(Edge *edge, VertexList *vertices, VertexList *poly_verts, size_t num_sides, unsigned padding);
bool _is_polygonslist_empty(PolyList polies);
void _add_poly(Poly *poly, PolyList *polygons);
void draw_polygons(PolyList polygons, VertexList *poly_verts, Color color);
void draw_poly(VertexList *poly_verts, Poly poly, Color color);
void free_poly(Poly *poly);
void _free_polygons(PolyList *polygons);
void poly_from_edge_list(PolyList *polygons, EdgeList *edges, VertexList *vertices, VertexList *poly_verts, unsigned padding);
void road_edges(PolyList *polygons, VertexList *vertices, EdgeList *road_edgs);
void _merge_poly(Poly *a, Poly *b, EdgeList *edges, VertexList *vertices);
bool _point_in_poly(Poly *poly, Vector2 point, VertexList list);
bool _point_inside_any_polygon(Vector2 point, EdgeList *edges, VertexList *list);
bool _edge_is_inside(Edge *edge, EdgeList *existing_edges, VertexList *list);
bool calculateEdgeIntersection(Edge *a, Edge *b, VertexList *list, Vector2 *result);
void _push_edge(EdgeList *edges, Edge *edge);
void _split_edge_at_all_intersections(Edge *edge, EdgeList *existing_edges, EdgeList *result, VertexList *list);
void _poly_push_merge(EdgeList *list, Poly *poly, VertexList *vertices);
EdgeList _get_poly_edges(Poly *poly); // Declaration for _get_poly_edges

int main(void)
{
    VertexList vertices = init_vertices(0);
    VertexList poly_vertices = init_vertices(0);

    EdgeList edges = _init_edges();
    EdgeList poly_egdes = _init_edges(); // This seems unused, consider removing
    EdgeList road = _init_edges();

    PolyList polygons = init_polygons();
    Poly tempPoly = {0};

    InitWindow(width, height, "figures");
    // glob vars
    Vector2 mousePos;
    int prev = -1;
    int curr = -1;
    static int sel = -1;
    bool was_dragging = false;

    SetTargetFPS(60);
    while (!WindowShouldClose())
    {
        mousePos = GetMousePosition();

        // add adj graph
        if (IsKeyPressed(KEY_ENTER))
        {
            prev = -1;
            curr = -1;
        }

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT))
        {
            sel = selectedVertex(vertices, mousePos, 7);
            if (sel >= 0)
            {
                prev = curr;
                curr = sel;
                vertices.vertices[sel].selected = true;
                was_dragging = false;
            }
            else
            {
                // Create a new vertex at the mouse position
                Vertex newVertex = {mousePos, true, true};
                addVertex(&vertices, newVertex);
                prev = curr;
                curr = vertices.count - 1;
            }

            if (prev >= 0)
            {
                _add_edge(&edges, &vertices, prev, curr, false);
                // Regenerate all polygons after adding edge
                _free_polygons(&polygons);
                polygons = init_polygons();
                free_vertices(&poly_vertices);
                poly_vertices = init_vertices(0);

                poly_from_edge_list(&polygons, &edges, &vertices, &poly_vertices, 40);

                // Rebuild the road network from all polygons
                free_edges(&road);
                road = _init_edges();
                Poly *current_poly_to_merge = polygons.start;
                while (current_poly_to_merge != NULL)
                {
                    _poly_push_merge(&road, current_poly_to_merge, &poly_vertices);
                    current_poly_to_merge = current_poly_to_merge->next;
                }
            }

            printf("( %d,%d )\n", curr, prev);
        }
        else if (IsMouseButtonDown(MOUSE_BUTTON_LEFT) && sel > -1)
        {
            // Update vertex position while dragging
            vertices.vertices[sel].loc = mousePos;
            was_dragging = true;
            tempPoly = poly_from_edge(edges.end, &vertices, &poly_vertices, 4, 50);
        }
        else if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT) && sel > -1 && was_dragging)
        {
            // Regenerate polygons and road only when finished dragging
            free_vertices(&poly_vertices);
            poly_vertices = init_vertices(0); // reset list
            _free_polygons(&polygons);
            polygons = init_polygons();
            poly_from_edge_list(&polygons, &edges, &vertices, &poly_vertices, 40);

            // Rebuild the road network from all polygons after dragging
            free_edges(&road);
            road = _init_edges();
            Poly *current_poly_to_merge = polygons.start;
            while (current_poly_to_merge != NULL)
            {
                _poly_push_merge(&road, current_poly_to_merge, &poly_vertices);
                current_poly_to_merge = current_poly_to_merge->next;
            }

            sel = -1;
            was_dragging = false;
            free_poly(&tempPoly);    // Free tempPoly's corners
            tempPoly.corners = NULL; // Ensure it's marked as freed
        }

        //------------- draw -------------//
        ClearBackground(RAYWHITE);
        BeginDrawing();

        drawVertices(vertices, 7);
        draw_edges(&edges, vertices);
        // draw_polygons(polygons, &poly_vertices, BLUE);
        draw_edges(&road, poly_vertices); // Draw the final road network

        if (was_dragging && edges.end != NULL && tempPoly.corners != NULL) // Check if tempPoly is valid
        {
            draw_poly(&poly_vertices, tempPoly, GRAY);
        }
        DrawFPS(10, 10);
        EndDrawing();
    }

    // Free memory when done
    _free_polygons(&polygons);
    free_edges(&edges);
    free_edges(&road); // Free road edges
    free_vertices(&vertices);
    free_vertices(&poly_vertices); // Free poly_vertices
    CloseWindow();

    return 0;
}

//---------------------------------------------------------
VertexList init_vertices(unsigned size)
{
    VertexList list;
    list.count = size;
    list.vertices = malloc(sizeof(Vertex) * size);
    if (list.vertices == NULL)
    {
        fprintf(stderr, "no space for a new list");
        exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < size; i++)
    {
        list.vertices[i] = (Vertex){
            (Vector2){(float)(rand() % width), (float)(rand() % height)}, // Cast to float
            false,
            true}; // Initialize editable
    }

    return list;
}

void drawVertices(VertexList list, int radius)
{
    for (size_t i = 0; i < list.count; i++)
    {
        drawVertex(list.vertices[i], radius);
    }
}

void drawVertex(Vertex v, int radius)
{
    if (v.selected)
    {
        DrawCircleV(v.loc, (float)radius + 3, GREEN); // Cast to float
    }
    DrawCircleV(v.loc, (float)radius, v.selected ? WHITE : GRAY); // Cast to float
    DrawCircleV(v.loc, (float)radius - 2, BLACK);                 // Cast to float
}

int selectedVertex(VertexList list, Vector2 pos, int radius)
{
    int index = -1;
    for (size_t i = 0; i < list.count; i++)
    {
        if (CheckCollisionCircles(list.vertices[i].loc, (float)radius, pos, 1)) // Cast to float
        {
            index = i;
        }
        list.vertices[i].selected = false; // Reset selection for all
    }
    return index;
}

void update_vertex(VertexList list, Vector2 newPos, int index)
{
    assert(index < list.count);
    list.vertices[index].loc = newPos;
}

void addVertex(VertexList *list, Vertex v)
{
    int old_count = list->count;
    list->count += 1;
    list->vertices = realloc(list->vertices, sizeof(Vertex) * (list->count));
    if (list->vertices == NULL)
    {
        fprintf(stderr, "Failed to reallocate memory for new vertex");
        exit(EXIT_FAILURE);
    }

    // Add the new vertex
    list->vertices[old_count] = v;
}

void pop_vertex(VertexList *list)
{
    if (list == NULL || list->count <= 0)
    {
        return;
    }

    list->count--;
}

void free_vertices(VertexList *list)
{
    if (list->vertices != NULL)
    {
        free(list->vertices);
        list->vertices = NULL;
    }
    list->count = 0;
}

//-----------------------------------
EdgeList _init_edges(void)
{
    EdgeList edges = (EdgeList){NULL, NULL, 0};
    return edges;
}

bool _is_empty_edgesList(EdgeList *edges)
{
    if (edges->begin == NULL && edges->end == NULL)
        return true;
    return false;
}

void _add_edge(EdgeList *edges, VertexList *list, int start, int end, bool overlap)
{
    Edge *edge = malloc(sizeof(Edge));
    if (edge == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for new edge\n");
        exit(EXIT_FAILURE);
    }
    edge->start = start;
    edge->end = end;
    edge->id = edges->size + 1;
    edge->next = NULL;

    if (start == end)
    {
        printf("Skip edge add: start and end vertices are the same.\n");
        free(edge);
        return;
    }

    // Check for duplicate edges
    Edge *current_check = edges->begin;
    while (current_check != NULL)
    {
        if ((current_check->start == start && current_check->end == end) ||
            (current_check->start == end && current_check->end == start))
        {
            printf("Skip edge add: duplicate edge (%d, %d).\n", start, end);
            free(edge);
            return;
        }
        current_check = current_check->next;
    }

    if (_is_empty_edgesList(edges))
    {
        edges->begin = edge;
        edges->end = edge;
        edges->size++;
    }
    else
    {
        Edge *current = edges->begin;
        Edge *prev = NULL;

        while (current != NULL)
        {
            // Check if end vertex lays on one of existing edges
            bool intercept = _is_vertex_on_edge(list->vertices[end], *current, *list);
            if (intercept)
            {
                // samething loading here
            }

            Vertex cross_point = {0};
            cross_point.editable = false;
            cross_point.selected = false;
            bool is_cross = _edges_cross(*current, *edge, *list, &cross_point);
            if (is_cross && !overlap)
            {
                // not needed
            }
            prev = current;
            current = current->next;
        }
        edges->end->next = edge;
        edges->end = edge;
        edges->size++;
    }
}

void free_edges(EdgeList *edges)
{
    if (edges == NULL)
    {
        return;
    }

    Edge *current = edges->begin;
    Edge *next;
    while (current != NULL)
    {
        next = current->next;
        free(current);
        current = next;
    };

    edges->begin = NULL;
    edges->end = NULL;
    edges->size = 0;
}

void draw_edges(EdgeList *edges, VertexList list)
{
    if (_is_empty_edgesList(edges))
    {
        return;
    }
    else
    {
        Edge *current = edges->begin;
        while (current)
        {
            DrawLineV(list.vertices[current->start].loc, list.vertices[current->end].loc, GRAY);
            current = current->next;
        }
    }
}

bool _is_vertex_on_edge(Vertex v, Edge edge, VertexList list)
{
    Vertex start_v = list.vertices[edge.start];
    Vertex end_v = list.vertices[edge.end];

    // Calculate distance from point to line segment
    float A = v.loc.x - start_v.loc.x;
    float B = v.loc.y - start_v.loc.y;
    float C = end_v.loc.x - start_v.loc.x;
    float D = end_v.loc.y - start_v.loc.y;

    float dot = A * C + B * D;
    float len_sq = C * C + D * D;
    float param = 0.0f;

    if (len_sq != 0) // Avoid division by zero
        param = dot / len_sq;

    float xx, yy;

    if (param < 0 || len_sq == 0)
    {
        xx = start_v.loc.x;
        yy = start_v.loc.y;
    }
    else if (param > 1) // After end
    {
        xx = end_v.loc.x;
        yy = end_v.loc.y;
    }
    else // On segment
    {
        xx = start_v.loc.x + param * C;
        yy = start_v.loc.y + param * D;
    }

    float dx = v.loc.x - xx;
    float dy = v.loc.y - yy;
    float distance = sqrt(dx * dx + dy * dy);

    return distance < 5.0f;
}

void _put_vertex_between(EdgeList *edges, Edge *edge, Edge *prev_edge, Vertex vertex, int vertex_index)
{
    // not really used tho;

    Edge *start_mid = malloc(sizeof(Edge));
    if (start_mid == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for start_mid edge\n");
        exit(EXIT_FAILURE);
    }
    start_mid->start = edge->start;
    start_mid->end = vertex_index;
    start_mid->id = edges->size + 1;
    start_mid->next = NULL;

    Edge *mid_end = malloc(sizeof(Edge));
    if (mid_end == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for mid_end edge\n");
        free(start_mid);
        exit(EXIT_FAILURE);
    }
    mid_end->start = vertex_index;
    mid_end->end = edge->end;
    mid_end->id = edges->size + 2;
    mid_end->next = edge->next;

    start_mid->next = mid_end;

    if (prev_edge == NULL)
    {
        edges->begin = start_mid;
    }
    else
    {
        prev_edge->next = start_mid;
    }

    if (edge == edges->end)
    {
        edges->end = mid_end;
    }
    free(edge);
    edges->size += 1;
}

bool _edges_cross(Edge e1, Edge e2, VertexList list, Vertex *result)
{
    Vector2 A = list.vertices[e1.start].loc;
    Vector2 B = list.vertices[e1.end].loc;
    Vector2 C = list.vertices[e2.start].loc;
    Vector2 D = list.vertices[e2.end].loc;

    const float tTop = (D.x - C.x) * (A.y - C.y) - (D.y - C.y) * (A.x - C.x);
    const float uTop = (C.y - A.y) * (A.x - B.x) - (C.x - A.x) * (A.y - B.y);
    const float bottom = (D.y - C.y) * (B.x - A.x) - (D.x - C.x) * (B.y - A.y);

    if (fabs(bottom) < 1e-6)
        return false;

    const float t = tTop / bottom;
    const float u = uTop / bottom;

    if (t > 1e-6 && t < 1.0f - 1e-6 && u > 1e-6 && u < 1.0f - 1e-6)
    {
        result->loc.x = A.x + t * (B.x - A.x);
        result->loc.y = A.y + t * (B.y - A.y);
        return true;
    }

    return false;
}

Split _split_edge(EdgeList *edges, Edge *edge, int vertex)
{
    // split edges

    Edge *first_part = malloc(sizeof(Edge));
    if (first_part == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for first_part\n");
        exit(EXIT_FAILURE);
    }
    first_part->start = edge->start;
    first_part->end = vertex;
    first_part->id = edge->id;
    first_part->next = NULL;

    Edge *second_part = malloc(sizeof(Edge));
    if (second_part == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for second_part\n");
        free(first_part);
        exit(EXIT_FAILURE);
    }
    second_part->start = vertex;
    second_part->end = edge->end;
    second_part->id = edge->id + 1;
    second_part->next = edge->next;

    first_part->next = second_part;

    return (Split){first_part, second_part};
}

// ----------------- polygons ------------- //

Poly init_poly(size_t shape)
{
    Poly pol = {0};
    pol.shape = shape;
    return pol;
}

Poly poly_from_edge(Edge *edge, VertexList *vertices, VertexList *poly_verts, size_t num_sides, unsigned padding)
{
    Poly poly = {0};
    poly.id = edge->id;
    poly.shape = num_sides;
    poly.corners = malloc(sizeof(int) * num_sides);
    if (poly.corners == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for poly corners\n");
        exit(EXIT_FAILURE);
    }

    Vector2 start = vertices->vertices[edge->start].loc;
    Vector2 end = vertices->vertices[edge->end].loc;

    Vector2 direction = {end.x - start.x, end.y - start.y};
    float length = Vector2Length(direction);

    if (length > 0)
    {
        direction = Vector2Normalize(direction);
    }

    Vector2 perp = {-direction.y * (float)padding / 2, direction.x * (float)padding / 2};

    // Store corner positions directly (no vertices added)
    Vertex a = {(Vector2){start.x - perp.x, start.y - perp.y}, false, false};
    Vertex b = {(Vector2){start.x + perp.x, start.y + perp.y}, false, false};
    Vertex c = {(Vector2){end.x + perp.x, end.y + perp.y}, false, false};
    Vertex d = {(Vector2){end.x - perp.x, end.y - perp.y}, false, false};

    // Add vertices to the poly_verts list and store their indices
    addVertex(poly_verts, a);
    poly.corners[0] = poly_verts->count - 1;
    addVertex(poly_verts, b);
    poly.corners[1] = poly_verts->count - 1;
    addVertex(poly_verts, c);
    poly.corners[2] = poly_verts->count - 1;
    addVertex(poly_verts, d);
    poly.corners[3] = poly_verts->count - 1;

    return poly;
}

void poly_from_edge_list(PolyList *polygons, EdgeList *edges, VertexList *vertices, VertexList *poly_verts, unsigned padding)
{
    if (_is_empty_edgesList(edges))
    {
        return;
    }

    Edge *current = edges->begin;
    while (current != NULL)
    {
        Poly *new_poly = malloc(sizeof(Poly));
        if (new_poly == NULL)
        {
            fprintf(stderr, "Failed to allocate memory for polygon\n");
            exit(EXIT_FAILURE);
        }

        *new_poly = poly_from_edge(current, vertices, poly_verts, 4, padding);
        new_poly->next = NULL;

        _add_poly(new_poly, polygons);
        current = current->next;
    }
}

void draw_poly(VertexList *poly_verts, Poly poly, Color color)
{
    for (size_t i = 0; i < poly.shape; i++)
    {
        int next = (int)((i + 1) % poly.shape); // (from 1 - shape)
        Vector2 start = poly_verts->vertices[poly.corners[i]].loc;
        Vector2 end = poly_verts->vertices[poly.corners[next]].loc;
        DrawLineEx(start, end, 3, color);
    }
}

bool _is_valid_poly(Poly *poly)
{
    if (poly == NULL || poly->corners == NULL || poly->shape == 0)
    {
        return false;
    }
    return true;
}

PolyList init_polygons(void)
{
    return (PolyList){NULL, NULL, 0};
}

bool _is_polygonslist_empty(PolyList polies)
{
    return (polies.start == NULL && polies.end == NULL && polies.count == 0);
}

void _add_poly(Poly *poly, PolyList *polygons)
{
    assert(poly != NULL);
    if (_is_polygonslist_empty(*polygons))
    {
        polygons->start = poly;
        polygons->end = poly;
    }
    else
    {
        polygons->end->next = poly;
        polygons->end = poly;
    }
    polygons->count++;
}

void draw_polygons(PolyList polygons, VertexList *poly_verts, Color color)
{
    if (_is_polygonslist_empty(polygons))
    {
        return;
    }

    Poly *current = polygons.start;
    while (current != NULL)
    {
        draw_poly(poly_verts, *current, color);
        current = current->next;
    }
}

EdgeList _get_poly_edges(Poly *poly)
{
    if (!_is_valid_poly(poly))
    {
        fprintf(stderr, "no valid polygon\n");
        assert(false);
    }

    EdgeList edges = _init_edges();

    if (poly->shape < 2)
    {
        return edges;
    }

    for (size_t i = 0; i < poly->shape; i++)
    {
        int start_corner = poly->corners[i];
        int end_corner = poly->corners[(i + 1) % poly->shape];

        Edge *edge = malloc(sizeof(Edge));
        if (edge == NULL)
        {
            fprintf(stderr, "no space for new edge\n");
            free_edges(&edges);
            exit(EXIT_FAILURE);
        }
        edge->start = start_corner;
        edge->end = end_corner;
        edge->id = poly->id * 1000 + (int)i;
        edge->next = NULL;

        _push_edge(&edges, edge);
    }

    return edges;
}

void free_poly(Poly *poly)
{
    if (poly != NULL && poly->corners != NULL)
    {
        free(poly->corners);
        poly->corners = NULL;
        poly->shape = 0;
        poly->id = 0;
    }
}

void _free_polygons(PolyList *polygons)
{
    if (_is_polygonslist_empty(*polygons))
    {
        return;
    }

    Poly *current = polygons->start;
    Poly *next = NULL;

    while (current != NULL)
    {
        next = current->next;
        free_poly(current);
        free(current);
        current = next;
    }

    polygons->start = NULL;
    polygons->end = NULL;
    polygons->count = 0;
}

void road_edges(PolyList *polygons, VertexList *vertices, EdgeList *road_edges)
{
    //
}

void _merge_poly(Poly *a, Poly *b, EdgeList *edges, VertexList *vertices)
{
    // This
}

bool _point_in_poly(Poly *poly, Vector2 point, VertexList list)
{
    if (!_is_valid_poly(poly))
        return false;

    int cross_count = 0;

    for (size_t i = 0; i < poly->shape; i++)
    {
        int next = (int)((i + 1) % poly->shape);

        Vector2 v1 = list.vertices[poly->corners[i]].loc;
        Vector2 v2 = list.vertices[poly->corners[next]].loc;
        // Raycast algorithm
        if (((v1.y > point.y) != (v2.y > point.y)) &&
            (point.x < (v2.x - v1.x) * (point.y - v1.y) / (v2.y - v1.y) + v1.x))
        {
            cross_count++;
        }
    }

    return (cross_count % 2) == 1; // Odd number of crossings = inside
}

// check if point inside *any* polygon represented by the given edges
bool _point_inside_any_polygon(Vector2 point, EdgeList *edges, VertexList *list)
{
    if (edges == NULL || _is_empty_edgesList(edges))
        return false;

    int crossings = 0;
    Edge *current = edges->begin;

    while (current != NULL)
    {
        Vector2 v1 = list->vertices[current->start].loc;
        Vector2 v2 = list->vertices[current->end].loc;

        if (((v1.y <= point.y && v2.y > point.y) || (v2.y <= point.y && v1.y > point.y)) &&
            (point.x < (v2.x - v1.x) * (point.y - v1.y) / (v2.y - v1.y) + v1.x))
        {
            crossings++;
        }
        current = current->next;
    }

    return (crossings % 2) == 1;
}

bool _edge_is_inside(Edge *edge, EdgeList *existing_edges, VertexList *list)
{
    if (edge == NULL || existing_edges == NULL)
        return false;

    // Check midpoint of the edge
    Vertex start = list->vertices[edge->start];
    Vertex end = list->vertices[edge->end];

    Vector2 midp = Vector2Scale(Vector2Add(start.loc, end.loc), 0.5f); //!!? Use 0.5f for midpoint

    return _point_inside_any_polygon(midp, existing_edges, list);
}

bool calculateEdgeIntersection(Edge *a, Edge *b, VertexList *list, Vector2 *result)
{
    if (a == NULL || b == NULL || result == NULL)
    {
        return false;
    }

    Vector2 A = list->vertices[a->start].loc;
    Vector2 B = list->vertices[a->end].loc;
    Vector2 C = list->vertices[b->start].loc;
    Vector2 D = list->vertices[b->end].loc;

    const float tTop = (D.x - C.x) * (A.y - C.y) - (D.y - C.y) * (A.x - C.x);
    const float uTop = (C.y - A.y) * (A.x - B.x) - (C.x - A.x) * (A.y - B.y);
    const float bottom = (D.y - C.y) * (B.x - A.x) - (D.x - C.x) * (B.y - A.y);

    if (fabs(bottom) < 1e-6)
        return false;

    const float t = tTop / bottom;
    const float u = uTop / bottom;

    if (t >= 0 && t <= 1 && u >= 0 && u <= 1)
    {
        result->x = A.x + t * (B.x - A.x);
        result->y = A.y + t * (B.y - A.y);

        float tolerance = 1.0f;
        if ((Vector2Distance(A, *result) < tolerance) ||
            (Vector2Distance(B, *result) < tolerance) ||
            (Vector2Distance(C, *result) < tolerance) ||
            (Vector2Distance(D, *result) < tolerance))
        {
            return false;
        }

        return true;
    }

    return false;
}

void _push_edge(EdgeList *edges, Edge *edge)
{
    if (edges == NULL || edge == NULL)
    {
        return;
    }

    edge->next = NULL;

    if (_is_empty_edgesList(edges))
    {
        edges->begin = edge;
        edges->end = edge;
        edges->size = 1;
    }
    else
    {
        edges->end->next = edge;
        edges->end = edge;
        edges->size++;
    }
}

void _split_edge_at_all_intersections(Edge *edge, EdgeList *existing_edges, EdgeList *result, VertexList *list)
{
    if (edge == NULL || existing_edges == NULL || result == NULL)
        return;

    typedef struct Intersection
    {
        int point_vertex_id;
        float t;
        struct Intersection *next;
    } Intersection;

    Intersection *intersections = NULL;
    int intersection_count = 0;

    Edge *current_existing_edge = existing_edges->begin;
    while (current_existing_edge != NULL)
    {
        if ((edge->start == current_existing_edge->start && edge->end == current_existing_edge->end) ||
            (edge->start == current_existing_edge->end && edge->end == current_existing_edge->start))
        {
            current_existing_edge = current_existing_edge->next;
            continue; // sk if same edge
        }
        if (edge->start == current_existing_edge->start || edge->start == current_existing_edge->end ||
            edge->end == current_existing_edge->start || edge->end == current_existing_edge->end)
        {
            current_existing_edge = current_existing_edge->next;
            continue;
        }

        Vector2 int_point_loc;
        if (calculateEdgeIntersection(edge, current_existing_edge, list, &int_point_loc))
        {
            Vertex intrPoint = {int_point_loc, false, false};
            addVertex(list, intrPoint);
            int intersection_vertex_id = list->count - 1;

            // Calculate parameter t for the current edge
            Vector2 edge_start_loc = list->vertices[edge->start].loc;
            Vector2 edge_end_loc = list->vertices[edge->end].loc;
            Vector2 dEdge = Vector2Subtract(edge_end_loc, edge_start_loc);
            float t_val = 0.0f;
            if (fabs(dEdge.x) > fabs(dEdge.y))
                t_val = (int_point_loc.x - edge_start_loc.x) / dEdge.x;
            else if (fabs(dEdge.y) > 1e-6) // Avoid division by zero for vertical lines
                t_val = (int_point_loc.y - edge_start_loc.y) / dEdge.y;
            else
                continue;
            if (t_val > 1e-6 && t_val < 1.0f - 1e-6)
            {
                // Check for duplicates (points very close to existing intersection points)
                bool is_duplicate = false;
                Intersection *check_dup = intersections;
                while (check_dup != NULL)
                {
                    // Check if point is already recorded or very close
                    if (Vector2Distance(list->vertices[check_dup->point_vertex_id].loc, int_point_loc) < 1.0f) // Tolerance for duplicates
                    {
                        is_duplicate = true;
                        // Use the existing vertex ID if it's a duplicate
                        intersection_vertex_id = check_dup->point_vertex_id;
                        break;
                    }
                    check_dup = check_dup->next;
                }

                if (!is_duplicate)
                {
                    Intersection *new_int = malloc(sizeof(Intersection));
                    if (new_int == NULL)
                    {
                        fprintf(stderr, "Failed to allocate memory for Intersection\n");
                        exit(EXIT_FAILURE);
                    }
                    new_int->point_vertex_id = intersection_vertex_id;
                    new_int->t = t_val;
                    new_int->next = NULL;

                    // Insert in sorted order by t value
                    if (intersections == NULL || t_val < intersections->t)
                    {
                        new_int->next = intersections;
                        intersections = new_int;
                    }
                    else
                    {
                        Intersection *curr_int_list = intersections;
                        while (curr_int_list->next != NULL && curr_int_list->next->t < t_val)
                            curr_int_list = curr_int_list->next;
                        new_int->next = curr_int_list->next;
                        curr_int_list->next = new_int;
                    }
                    intersection_count++;
                }
            }
        }
        current_existing_edge = current_existing_edge->next;
    }

    // Create edge segments based on collected intersections
    int current_segment_start_id = edge->start;
    Intersection *curr_int = intersections;
    int segment_idx = 0;

    while (curr_int != NULL)
    {
        // Create segment from current_segment_start_id to intersection point
        Edge *segment = malloc(sizeof(Edge));
        if (segment == NULL)
        {
            fprintf(stderr, "Failed to allocate memory for segment\n");
            exit(EXIT_FAILURE);
        }
        segment->start = current_segment_start_id;
        segment->end = curr_int->point_vertex_id;
        segment->id = edge->id * 1000 + segment_idx++;
        segment->next = NULL;

        if (!_edge_is_inside(segment, existing_edges, list))
        {
            _push_edge(result, segment);
        }
        else
        {
            free(segment);
        }

        current_segment_start_id = curr_int->point_vertex_id;
        curr_int = curr_int->next;
    }

    Edge *final_segment = malloc(sizeof(Edge));
    if (final_segment == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for final_segment\n");
        exit(EXIT_FAILURE);
    }
    final_segment->start = current_segment_start_id;
    final_segment->end = edge->end;
    final_segment->id = edge->id * 1000 + segment_idx;
    final_segment->next = NULL;

    if (!_edge_is_inside(final_segment, existing_edges, list))
    {
        _push_edge(result, final_segment);
    }
    else
    {
        free(final_segment);
    }

    while (intersections != NULL)
    {
        Intersection *temp = intersections;
        intersections = intersections->next;
        free(temp);
    }
}

void _poly_push_merge(EdgeList *road_list, Poly *new_poly, VertexList *vertices)
{
    if (new_poly == NULL || road_list == NULL)
    {
        fprintf(stderr, "empty list or polygon\n");
        return;
    }

    EdgeList new_poly_edges = _get_poly_edges(new_poly);

    if (_is_empty_edgesList(road_list))
    {
        road_list->begin = new_poly_edges.begin;
        road_list->end = new_poly_edges.end;
        road_list->size = new_poly_edges.size;
        new_poly_edges.begin = NULL;
        new_poly_edges.end = NULL;
        new_poly_edges.size = 0;
        return;
    }

    EdgeList new_poly_clipped_edges = _init_edges();
    EdgeList existing_road_filtered_edges = _init_edges();
    Edge *current_new_poly_edge = new_poly_edges.begin;
    while (current_new_poly_edge != NULL)
    {
        _split_edge_at_all_intersections(current_new_poly_edge, road_list, &new_poly_clipped_edges, vertices);
        current_new_poly_edge = current_new_poly_edge->next;
    }

    Edge *current_existing_road_edge = road_list->begin;
    while (current_existing_road_edge != NULL)
    {
        _split_edge_at_all_intersections(current_existing_road_edge, &new_poly_edges, &existing_road_filtered_edges, vertices);
        current_existing_road_edge = current_existing_road_edge->next;
    }

    EdgeList final_road_edges = _init_edges();

    Edge *edge_to_add = existing_road_filtered_edges.begin;
    while (edge_to_add != NULL)
    {
        Edge *next_edge = edge_to_add->next;
        edge_to_add->next = NULL;
        _push_edge(&final_road_edges, edge_to_add);
        edge_to_add = next_edge;
    }

    // Add new polygon's clipped edges
    edge_to_add = new_poly_clipped_edges.begin;
    while (edge_to_add != NULL)
    {
        Edge *next_edge = edge_to_add->next;
        edge_to_add->next = NULL; // Detach from old list
        _push_edge(&final_road_edges, edge_to_add);
        edge_to_add = next_edge;
    }

    free_edges(road_list);
    road_list->begin = final_road_edges.begin;
    road_list->end = final_road_edges.end;
    road_list->size = final_road_edges.size;

    new_poly_edges.begin = NULL;
    new_poly_edges.end = NULL;
    new_poly_edges.size = 0;
    free_edges(&new_poly_edges);

    existing_road_filtered_edges.begin = NULL;
    existing_road_filtered_edges.end = NULL;
    existing_road_filtered_edges.size = 0;
    free_edges(&existing_road_filtered_edges);

    new_poly_clipped_edges.begin = NULL;
    new_poly_clipped_edges.end = NULL;
    new_poly_clipped_edges.size = 0;
    free_edges(&new_poly_clipped_edges);
}
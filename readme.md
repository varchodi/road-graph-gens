# Road Network Generator

This C program, built with the Raylib library, allows users to interactively create a road-like network by defining polygons. The system automatically merges overlapping polygon boundaries to form a continuous "road" or pathway.

## Features

- **Interactive Vertex and Edge Creation:** Click to add vertices and connect them with edges to form a base graph.
- **Polygon Generation:** Polygons are automatically generated from the defined edges, representing road segments.
- **Dynamic Road Merging:** Overlapping or adjacent polygons are merged to create a unified "road" network, effectively performing a union operation on the polygon boundaries.
- **Real-time Visualization:** See your road network evolve in real-time with Raylib's drawing capabilities.
- **Vertex Dragging:** Click and drag existing vertices to modify the underlying graph and observe the road network update dynamically.

## How it Works

The core of this system lies in the `_poly_push_merge` function. When a new polygon is created (or existing ones are modified), this function performs the following steps:

1.  **Extract Polygon Edges:** It first extracts the individual edges that make up the new polygon.
2.  **Split and Filter New Edges:** Each edge of the _new_ polygon is checked against the _existing_ road network. Only the portions of the new polygon's edges that lie _outside_ the current road network are kept.
3.  **Split and Filter Existing Road Edges:** Conversely, each edge of the _existing_ road network is checked against the _new_ polygon. Only the portions of the existing road edges that lie _outside_ the new polygon are retained.
4.  **Combine and Replace:** The filtered segments from both the new polygon and the existing road network are combined to form the updated, merged road network. This process effectively creates a union of the polygon boundaries.

This iterative merging process ensures that as you add more polygons, the `road` `EdgeList` always represents the outer boundary of all combined polygon shapes.

## Getting Started

### Prerequisites

- **Raylib Library:** This project heavily relies on Raylib for graphics and input. You'll need to have Raylib installed and configured for your development environment.
  - [Raylib Installation Guide](https://github.com/raysan5/raylib/wiki/Setup-and-Installation)

### Building the Project

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/varchodi/road-graph-gens.git](https://github.com/varchodi/road-graph-gens.git)
    cd road-network-generator
    ```
2.  **Compile with Raylib:**
    Assuming you have Raylib set up, you can compile the `main.c` file using a C compiler (like GCC) and linking with Raylib.

    On Linux/macOS:

    ```bash
    gcc main.c -o road_shit -lraylib -lm
    ```

    On Windows (using MinGW/MSYS2):

    ```bash
    gcc main.c -o road_shit.exe -lraylib -lopengl32 -lgdi32 -lwinmm -static
    ```

    _(Note: The exact linking flags might vary depending on your Raylib installation.)_

### Running the Application

After successful compilation, run the executable:

```bash
./road_shit
```

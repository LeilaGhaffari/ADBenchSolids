import matplotlib.pyplot as plt
import networkx as nx

def plot_ast_graph():
    ast_graph = nx.DiGraph()
    ast_graph.add_edges_from([
        ("x", "pow2"),
        ("pow2", "*"),
        ("y", "*"),
        ("x", "sin"),
        ("sin", "+"),
        ("*", "+"),
        ("+", "f"),
    ])

    pos = {
        "f": (1, 10),
        "+": (1, 5),
        "sin": (.5, 2),
        "*": (1.5, 2),
        "pow2": (1, -.5),
        "x": (0.25, -1),
        "y": (2, -1),
    }

    node_labels = {
        "x": "x",
        "y": "y",
        "pow2": "z₁ = x²",
        "*": "z₂ = z₁·y",
        "sin": "z₃ = sin(x)",
        "+": "z₄ = z₂ + z₃",
        "f": "f = z₄",
    }

    node_colors = []
    for node in ast_graph.nodes:
        if node == "f":
            node_colors.append("lightgreen")
        elif node == "x" or node == "y":
            node_colors.append("orange")
        else:
            node_colors.append("skyblue")

    nx.draw(
        ast_graph, pos, with_labels=True, node_color=node_colors, edge_color='black',
        node_size=3000, font_size=10, font_weight='bold', arrowsize=20
    )

    nx.draw_networkx_edge_labels(
        ast_graph, pos,
        edge_labels={
            ("pow2", "*"): "z₁ = x²",
            ("sin", "+"): "z₃ = sin(x)",
            ("*", "+"): "z₂ = z₁·y",
            ("+", "f"): "f = z₂ + z₃",
        },
        font_color='red',
        label_pos=0.5,
        rotate=False,
        font_size=10
    )



    plt.title("Abstract Syntax Tree (AST)", fontsize=14)
    plt.tight_layout()
    plt.savefig("ast_graph.png")
    plt.close()


def plot_computational_graph():
    forward_graph = nx.DiGraph()
    forward_graph.add_edges_from([
        ("x", "pow2"),
        ("pow2", "*"),
        ("y", "*"),
        ("x", "sin"),
        ("sin", "+"),
        ("*", "+"),
        ("+", "f"),
    ])

    pos = {
        "f": (1, 10),
        "+": (1, 5),
        "sin": (.5, 2),
        "*": (1.5, 2),
        "pow2": (1, -.5),
        "x": (0.25, -1),
        "y": (2, -1),
    }

    node_labels = {
        "x": "x",
        "y": "y",
        "pow2": "z₁ = x²",
        "*": "z₂ = z₁·y",
        "sin": "z₃ = sin(x)",
        "+": "z₄ = z₂ + z₃",
        "f": "f = z₄",
    }

    node_colors = []
    for node in forward_graph.nodes:
        if node == "f":
            node_colors.append("lightgreen")
        elif node == "x" or node == "y":
            node_colors.append("orange")
        else:
            node_colors.append("skyblue")

    nx.draw(
        forward_graph, pos, node_color=node_colors, edge_color='black',
        node_size=3000, font_size=10, font_weight='bold', arrowsize=20
    )
    nx.draw_networkx_labels(forward_graph, pos, labels=node_labels, font_size=10, font_weight='bold')

    nx.draw_networkx_edge_labels(
        forward_graph, pos,
        edge_labels={
            ("x", "pow2"): "$dx$",
            ("pow2", "*"): "$dz_1=2xdx$",
            ("y", "*"): "$dy$",
            ("x", "sin"): "$dx$",
            ("sin", "+"): "$dz_3=cos(x)dx$",
            ("*", "+"): "$dz_2=dz_1y+x^2 dy$",
            ("+", "f"): "$df=dz_2+dz_3$",
        },
        font_color='red',
        label_pos=0.5,
        rotate=False
    )

    plt.title("Forward Mode Computational Graph (Flipped)", fontsize=14)
    plt.tight_layout()
    plt.savefig("forward_mode_graph.png")
    plt.close()


    # Reverse Mode Graph
    reverse_graph = nx.DiGraph()
    reverse_graph.add_edges_from([
        ("$\\bar{f}=1$", "sin"),
        ("$\\bar{f}=1$", "*"),
        ("sin", "+"),
        ("*", "pow2"),
        ("*", "$\\bar{y}$"),
        ("pow2", "+"),
        ("+", "$\\bar{x}$"),
    ])
    pos = {
        "$\\bar{f}=1$": (1, 10),
        "sin": (.5, 6),
        "*": (1.5, 6),
        "pow2": (1.1, 2.5),
        "+": (.25, 1),
        "$\\bar{x}$": (-.25, -1),
        "$\\bar{y}$": (2, -1),
    }

    # Define node colors for Reverse Mode
    node_colors_reverse = []
    for node in reverse_graph.nodes:
        if node == "$\\bar{f}=1$":
            node_colors_reverse.append("orange")
        elif node == "$\\bar{x}$" or node == "$\\bar{y}$":
            node_colors_reverse.append("lightgreen")
        else:
            node_colors_reverse.append("skyblue")

    nx.draw(
        reverse_graph, pos, with_labels=True, node_color=node_colors_reverse, edge_color='black',
        node_size=3000, font_size=10, font_weight='bold', arrowsize=20
    )

    nx.draw_networkx_edge_labels(
        reverse_graph, pos,
        edge_labels={
            ("$\\bar{f}=1$", "sin"): "$\\bar{z}_3 = \\bar{f} \cdot \\frac{\partial f}{\partial z_3}=1$",
            ("$\\bar{f}=1$", "*"): "$\\bar{z}_2 = \\bar{f} \cdot \\frac{\partial f}{\partial z_2}=1$",
            ("sin", "+"): "$\\bar{z}_3 \cdot \\frac{\partial z_3}{\partial x}=cos(x)$",
            ("*", "pow2"): "$\\bar{z}_1 = \\bar{z}_2 \cdot \\frac{\partial z_2}{\partial z_1}=y$",
            ("*", "$\\bar{y}$"): "$\\bar{y} = \\bar{z}_2 \cdot \\frac{\partial z_2}{\partial y} = z_1 = x^2$",
            ("pow2", "+"): "$\\bar{z}_1 \cdot \\frac{\partial z_1}{\partial x}=2xy$",
        },
        font_color='red',
        label_pos=0.5,
        rotate=False
    )

    plt.title("Reverse Mode Computational Graph (Flipped)", fontsize=14)
    plt.tight_layout()
    plt.savefig("reverse_mode_graph.png")
    plt.close()


if __name__ == "__main__":
    plot_computational_graph()
    plot_ast_graph()
    print("Computational graphs saved as 'forward_mode_graph.png', 'reverse_mode_graph.png', and 'ast_graph.png'.")


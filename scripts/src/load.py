import yaml
import os

# Get the repository root (2 levels up from this file)
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _resolve_paths(node):
    """Resolve relative paths by prefixing with REPO_ROOT; keep absolute paths as-is."""
    if isinstance(node, dict):
        return {k: _resolve_paths(v) for k, v in node.items()}
    if isinstance(node, list):
        return [_resolve_paths(v) for v in node]
    if isinstance(node, str):
        return node if os.path.isabs(node) else os.path.normpath(os.path.join(REPO_ROOT, node))
    return node


# Load config
def load_config():
    config_path = os.path.join(REPO_ROOT, "config.yaml")
    with open(config_path) as f:
        config = yaml.safe_load(f)
    return _resolve_paths(config)


def get_root():
    return REPO_ROOT


if __name__ == "__main__":
    config = load_config()
    print(config)
# pfl_io: scene save/load (JSON), the built-in example scenes, CSV data export and the GIF encoder.
# No graphics, so the unit tests link it too.
pfl_embed_resources(pfl_examples
    BASE_DIR "${PROJECT_SOURCE_DIR}"
    NAMESPACE examples
    FILES
        examples/uniform_45deg.json
        examples/rankine_half_body.json
        examples/rankine_oval.json
        examples/cylinder.json
        examples/lifting_cylinder.json
        examples/source_near_wall.json
        examples/channel_source.json
        examples/source_sink_vortex.json
        examples/corner_flow_90deg.json
        examples/flat_plate_edge.json)

add_library(pfl_io STATIC
    scene_io.hpp scene_io.cpp
    examples.hpp examples.cpp
    data_export.hpp data_export.cpp
    verification_checks.hpp verification_checks.cpp
    gif_writer.hpp gif_writer.cpp)
target_include_directories(pfl_io PUBLIC "${PROJECT_SOURCE_DIR}/src")
target_link_libraries(pfl_io PUBLIC pfl_core PRIVATE nlohmann_json::nlohmann_json pfl_examples)
pfl_set_warnings(pfl_io)

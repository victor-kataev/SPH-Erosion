#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"


#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <math.h>
#include <iostream>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


#include "grid.h"
#include "shader.h"
#include "camera.h"
#include "fluid_system.h"
#include "mesh.h"
#include "shapes.h"

#define SCREEN_WIDTH 1440
#define SCREEN_HEIGHT 900
#define GL_PI 3.1415f

GLint w = SCREEN_WIDTH;
GLint h = SCREEN_HEIGHT;
float deltaTime = 0;
float lastFrame = 0;

Camera camera(glm::vec3(7.0, 256.0, 10.0));
float lastX = SCREEN_WIDTH / 2.0f;
float lastY = SCREEN_HEIGHT / 2.0f;
bool firstMouse = true;
bool disabled = true;
bool pressed_before = false;


void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
void checkGLError(const char* where, int line);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);

#define CHECK_GL_ERROR() do { checkGLError(__FUNCTION__, __LINE__); } while (0)

FluidSystemSPH fluidsph;


int main()
{
    if (!glfwInit())
        return -1;

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Erosion", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create glfw window\n";
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);


    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    bool err = glewInit() != GLEW_OK;
    if (err)
    {
        fprintf(stderr, "Failed to initialize OpenGL loader!\n");
        return 1;
    }

    Shader shader("vertex.glsl", "fragment.glsl");
    //Shader surface_shader("vertex.glsl", "fragment_surface.glsl");
  
    char picture_path[100];
    strcpy_s(picture_path, "pumba_gray.png");

    int width, height, channels;
    unsigned char* img = stbi_load(picture_path, &width, &height, &channels, 1);
    if (img == NULL)
    {
        printf("Error in loading the image\n");
        exit(1);
    }
    printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n", width, height, channels);
    
    float dimensions[3] = { 50, 355, 50};
    Grid grid(dimensions[0], dimensions[1], dimensions[2]);
    grid.LoadHeightfield(img);
    stbi_image_free(img);

    grid.UpdateGrid((int)dimensions[0], (int)dimensions[1], (int)dimensions[2]);

    
    std::vector<float> surface_verts = grid.GetSurfaceParts();
    //size_t surface_size = grid.GetSurfacePartsSize();
    std::vector<unsigned int> indices = grid.GetIndices();
    //size_t indices_size = grid.GetIndicesSize();
    //std::vector<float> fluid_verts = grid.GetFluidParts();
    //size_t fluid_size = grid.GetFluidPartsSize();
    

    float vertices[] = {
        1, 0, 1,
        1, 0, 2,
        -1, 0, 1,
        -1, 0, 2
    };

    unsigned int surfaceVAO, surfaceVBO, surfaceEBO;
    unsigned int fluidVAO, fluidVBO;

    glGenVertexArrays(1, &surfaceVAO);
    glGenBuffers(1, &surfaceVBO);
    glGenBuffers(1, &surfaceEBO);
    glBindVertexArray(surfaceVAO);
    glBindBuffer(GL_ARRAY_BUFFER, surfaceVBO);
    glBufferData(GL_ARRAY_BUFFER, surface_verts.size() * sizeof(float), &surface_verts[0], GL_STATIC_DRAW);
    //glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, surfaceEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);

    
   /* glGenVertexArrays(1, &fluidVAO);
    glGenBuffers(1, &fluidVBO);
    glBindVertexArray(fluidVAO);
    glBindBuffer(GL_ARRAY_BUFFER, fluidVBO);
    glBufferData(GL_ARRAY_BUFFER, fluid_verts.size() * sizeof(float), &fluid_verts[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glEnableVertexAttribArray(0);*/
    

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    
    // Setup Dear ImGui style
    ImGui::StyleColorsDark();

    // Setup Platform/Renderer backends
    const char* glsl_version = "#version 330";
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);

   
    bool show_demo_window = true;
    bool show_another_window = false;
    ImVec4 clear_color = ImVec4(0.6f, 0.0f, 0.0f, 1.00f);

    //fluidsph.SetOrigin(glm::vec3(32.9, 124.5, 43.2));
    fluidsph.SetOrigin(glm::vec3(0.0));
    //fluidsph.SetOrigin(glm::vec3(30.0, 255.5, 30.2));
    //camera.PlaceTo(glm::vec3(33.0, 125.0, 43.5));
    camera.PlaceTo(glm::vec3(0.0, 0.0, 2.0));
    fluidsph.Initialize(2);
    Shape shape;
    shape.CreateCube();
    shape.CreateBowl();

    Hemisphere hsphere(30, 30, 1, glm::vec3(0.0, 0.0, 0.0));
    Mesh terrain(surface_verts, indices);
    Mesh mHsphere(hsphere.GetVertices(), hsphere.GetIndices());



    glEnable(GL_DEPTH_TEST);

    //render loop
    while (!glfwWindowShouldClose(window))
    {
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        processInput(window);
        glClearColor(0.3, 0.3, 0.3, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        
        
        if (!disabled)
        {
            ImGui::Begin("New window");

            ImGui::ColorEdit3("clear color", (float*)&clear_color);
            ImGui::InputFloat3("dimensions", dimensions);
            if (ImGui::Button("Update grid"))
            {
                grid.UpdateGrid((int)dimensions[0], (int)dimensions[1], (int)dimensions[2]);
                surface_verts = grid.GetSurfaceParts();
                //surface_size = grid.GetSurfacePartsSize();
                indices = grid.GetIndices();
                //indices_size = grid.GetIndicesSize();

                glDeleteBuffers(1, &surfaceVBO);
                glDeleteBuffers(1, &surfaceEBO);

                glBindVertexArray(surfaceVAO);
                glGenBuffers(1, &surfaceVBO);
                glGenBuffers(1, &surfaceEBO);
                glBindBuffer(GL_ARRAY_BUFFER, surfaceVBO);
                glBufferData(GL_ARRAY_BUFFER, surface_verts.size() * sizeof(float), &surface_verts[0], GL_STATIC_DRAW);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, surfaceEBO);
                glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);
                glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
                glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
                glEnableVertexAttribArray(0);
            }
            ImGui::InputText("picture path", picture_path, 100);
            if (ImGui::Button("Load picture"))
            {
                img = stbi_load(picture_path, &width, &height, &channels, 1);
                if (img == NULL)
                {
                    printf("Error in loading the image\n");
                    exit(1);
                }
                printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n", width, height, channels);

                grid.LoadHeightfield(img);
                stbi_image_free(img);

                grid.UpdateGrid((int)dimensions[0], (int)dimensions[1], (int)dimensions[2]);
                surface_verts = grid.GetSurfaceParts();
                //surface_size = grid.GetSurfacePartsSize();
                indices = grid.GetIndices();
                //indices_size = grid.GetIndicesSize();

                glDeleteBuffers(1, &surfaceVBO);
                glDeleteBuffers(1, &surfaceEBO);

                glBindVertexArray(surfaceVAO);
                glGenBuffers(1, &surfaceVBO);
                glGenBuffers(1, &surfaceEBO);
                glBindBuffer(GL_ARRAY_BUFFER, surfaceVBO);
                glBufferData(GL_ARRAY_BUFFER, surface_verts.size() * sizeof(float), &surface_verts[0], GL_STATIC_DRAW);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, surfaceEBO);
                glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);
                glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
                glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
                glEnableVertexAttribArray(0);
            }

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
            ImGui::End();
        }

        
        
        fluidsph.Run(grid);//<------------------------

        glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)SCREEN_WIDTH / SCREEN_HEIGHT, 0.01f, 1000.0f);
        glm::mat4 view = camera.GetViewMatrix();
        glm::mat4 model = glm::mat4(1.0f);
        

        shader.use();
        shader.setMat4("projection", projection);
        shader.setMat4("view", view);
        shader.setMat4("model", model);
        shader.setVec3("dirLight.dir", glm::vec3(-0.1, -0.8, -0.8));
        shader.setVec3("dirLight.color", glm::vec3(1.0));
        shader.setVec3("viewerPos", camera.Position);
        shader.setVec3("material.ka", glm::vec3(0.2f));
        shader.setVec3("material.kd", glm::vec3(0.7f));
        shader.setVec3("material.ks", glm::vec3(1.0));
        shader.setFloat("material.ksh", 4.0f);


        shader.setVec3("myColor", glm::vec3(clear_color.x, clear_color.y, clear_color.z));
        glBindVertexArray(surfaceVAO);
        glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);


        //shader.setVec3("myColor", glm::vec3(1.0, 1.0, 0.0));
        //shape.DrawCube();

        shader.setVec3("myColor", glm::vec3(0.7, 0.8, 0.5));
        shape.DrawBowl();

        //shader.setVec3("myColor", glm::vec3(0.0, 1.0, 1.0));
        //hsphere.Draw();


        shader.setVec3("myColor", glm::vec3(0.0, 0.0, 1.0));
        fluidsph.Draw(shader);
        

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());


        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    //glDeleteVertexArrays(1, &surfaceVAO);
    //glDeleteBuffers(1, &surfaceVBO);
    
    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    
    glfwTerminate();


    return 0;
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}


bool pause = false;
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
            camera.ProcessKeyboard(FORWARD, deltaTime * 10);
        else if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
            camera.ProcessKeyboard(FORWARD, deltaTime * 0.1);
        else
            camera.ProcessKeyboard(FORWARD, deltaTime);

    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
            camera.ProcessKeyboard(BACKWARD, deltaTime * 10);
        else if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
            camera.ProcessKeyboard(BACKWARD, deltaTime * 0.1);
        else
            camera.ProcessKeyboard(BACKWARD, deltaTime);

    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
            camera.ProcessKeyboard(LEFT, deltaTime * 10);
        else if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
            camera.ProcessKeyboard(LEFT, deltaTime * 0.1);
        else
            camera.ProcessKeyboard(LEFT, deltaTime);

    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
            camera.ProcessKeyboard(RIGHT, deltaTime * 10);
        else if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
            camera.ProcessKeyboard(RIGHT, deltaTime * 0.1);
        else
            camera.ProcessKeyboard(RIGHT, deltaTime);

    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
            camera.ProcessKeyboard(UP, deltaTime * 10);
        else if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
            camera.ProcessKeyboard(UP, deltaTime * 0.1);
        else
            camera.ProcessKeyboard(UP, deltaTime);

    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
            camera.ProcessKeyboard(DOWN, deltaTime * 10);
        else if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
            camera.ProcessKeyboard(DOWN, deltaTime * 0.1);
        else
            camera.ProcessKeyboard(DOWN, deltaTime);

    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
        pressed_before = true;
    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_RELEASE && pressed_before)
    {
        pressed_before = false;
        if (disabled)
        {
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
            glfwSetCursorPosCallback(window, NULL);
            disabled = false;
            firstMouse = true;
        }
        else
        {
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
            glfwSetCursorPosCallback(window, mouse_callback);
            disabled = true;
        }
    }

    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS)
        pause = true;
    if (glfwGetKey(window, GLFW_KEY_F) == GLFW_RELEASE && pause)
    {
        if (fluidsph.GetDeltaTime() == 0)
            fluidsph.SetDeltaTime(0.01f);
        else
            fluidsph.SetDeltaTime(0.0f);
        pause = false;
    }

    if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS)
        fluidsph.PrintCoords();
    
    /*std::cout << "Position: "
        << camera.GetPosition().x << ' ' << camera.GetPosition().y << ' ' << camera.GetPosition().z
        << std::endl
        << "Front: "
        << camera.GetFront().x << ' ' << camera.GetFront().y << ' ' << camera.GetFront().z << std::endl;*/
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
    {
        fluidsph.Reset();
        std::cout << "Particles reseted\n";
    }

    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        fluidsph.AddParticles(125);
        std::cout << "Added 125 particles\n";
    }
}

void checkGLError(const char* where, int line)
{
    GLenum err = glGetError();
    if (err == GL_NONE)
        return;

    std::string errString = "<unknown>";
    switch (err) {
    case GL_INVALID_ENUM:
        errString = "GL_INVALID_ENUM";
        break;
    case GL_INVALID_VALUE:
        errString = "GL_INVALID_VALUE";
        break;
    case GL_INVALID_OPERATION:
        errString = "GL_INVALID_OPERATION";
        break;
    case GL_INVALID_FRAMEBUFFER_OPERATION:
        errString = "GL_INVALID_FRAMEBUFFER_OPERATION";
        break;
    case GL_OUT_OF_MEMORY:
        errString = "GL_OUT_OF_MEMORY";
        break;
    default:;
    }
    if (where == 0 || *where == 0)
        std::cerr << "GL error occurred: " << errString << std::endl;
    else
        std::cerr << "GL error occurred in " << where << ":" << line << ":" << errString << std::endl;
}

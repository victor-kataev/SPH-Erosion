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

#define SCREEN_WIDTH 1024
#define SCREEN_HEIGHT 768
#define GL_PI 3.1415f

GLint w = SCREEN_WIDTH;
GLint h = SCREEN_HEIGHT;
float deltaTime = 0;
float lastFrame = 0;

Camera camera(glm::vec3(100.0, 100.0, 103.0));
float lastX = SCREEN_WIDTH / 2.0f;
float lastY = SCREEN_HEIGHT / 2.0f;
bool firstMouse = true;
bool disabled = true;
bool pressed_before = false;


void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void framebuffer_size_callback(GLFWwindow* window, int width, int height);

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
        
    /*std::cout << "Position: "
        << camera.GetPosition().x << ' ' << camera.GetPosition().y << ' ' << camera.GetPosition().z
        << std::endl
        << "Front: "
        << camera.GetFront().x << ' ' << camera.GetFront().y << ' ' << camera.GetFront().z << std::endl;*/
}

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
    

    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    bool err = glewInit() != GLEW_OK;
    if (err)
    {
        fprintf(stderr, "Failed to initialize OpenGL loader!\n");
        return 1;
    }

    Shader shader("vertex.glsl", "fragment.glsl");
  
    char picture_path[100];
    strcpy_s(picture_path, "lena_gray.png");

    int width, height, channels;
    unsigned char* img = stbi_load(picture_path, &width, &height, &channels, 1);
    if (img == NULL)
    {
        printf("Error in loading the image\n");
        exit(1);
    }
    printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n", width, height, channels);
    
    float dimensions[3] = { 1, 1, 1 };
    Grid grid;
    grid.LoadHeightfield(img);
    stbi_image_free(img);

    grid.UpdateGrid((int)dimensions[0], (int)dimensions[1], (int)dimensions[2]);
    
    float* verts = grid.GetCoords();
    size_t size = grid.GetCoordsSize();


    unsigned int VBO, VAO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, size * sizeof(float), verts, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glEnableVertexAttribArray(0);


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
    ImVec4 clear_color = ImVec4(1.0f, 0.0f, 0.0f, 1.00f);

    while (!glfwWindowShouldClose(window))
    {
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        processInput(window);
        glClearColor(0.2, 0.3, 0.4, 1.0);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        

        //imguiSetUp();

        shader.use();

        shader.setVec3("myColor", glm::vec3(clear_color.x, clear_color.y, clear_color.z));

        glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)SCREEN_WIDTH / SCREEN_HEIGHT, 0.1f, 1000.0f);
        shader.setMat4("projection", projection);
        glm::mat4 view = camera.GetViewMatrix();
        shader.setMat4("view", view);

        glm::mat4 model = glm::mat4(1.0f);
        shader.setMat4("model", model);

        glBindVertexArray(VAO);

       // ImGui::ShowDemoWindow(&show_demo_window);


        if (!disabled)
        {
            ImGui::Begin("New window");                          // Create a window called "Hello, world!" and append into it.

            ImGui::ColorEdit3("clear color", (float*)&clear_color); // Edit 3 floats representing a color
            ImGui::InputFloat3("dimensions", dimensions);
            if (ImGui::Button("Update grid"))
            {
                grid.UpdateGrid((int)dimensions[0], (int)dimensions[1], (int)dimensions[2]);
                verts = grid.GetCoords();
                size = grid.GetCoordsSize();
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
                verts = grid.GetCoords();
                size = grid.GetCoordsSize();
            }

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
            ImGui::End();

        }

        glBufferData(GL_ARRAY_BUFFER, size * sizeof(float), verts, GL_STATIC_DRAW);
        glDrawArrays(GL_POINTS, 0, size);

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    
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
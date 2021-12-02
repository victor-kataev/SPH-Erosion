#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"


#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <math.h>
#include <iostream>
#include <stdio.h>
#include <Windows.h>
#include <sstream>


#include "grid.h"
#include "camera.h"
#include "fluid_system.h"
#include "mesh.h"
#include "shapes.h"
//#include "GLToMovie.h"
//#include <Vfw.h>

#define SCREEN_WIDTH 1440
#define SCREEN_HEIGHT 900
#define GL_PI 3.1415f

GLint w = SCREEN_WIDTH;
GLint h = SCREEN_HEIGHT;
float deltaTime = 0;
float lastFrame = 0;
int g_part_id = 0;


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
void pixelsToBmp(const char* filename, const unsigned char* pixels);

GLFWwindow* init();
void UIinit(GLFWwindow* window);
void UIshutdown();
void UIbegin();
void UIend();



#define CHECK_GL_ERROR() do { checkGLError(__FUNCTION__, __LINE__); } while (0)

FluidSystemSPH fluidsph;


int main()
{
    GLFWwindow* window = init();

    Shader shader("vertex.glsl", "fragment.glsl");
  
    char picture_path[100];
    //strcpy_s(picture_path, "lena_gray.png");
    strcpy_s(picture_path, "pumba_gray.png");
    glm::vec3 dimensions = { 50, 255, 50 };
    Grid grid(picture_path, dimensions);
        
    ImVec4 clear_color = grid.GetColor();

    fluidsph.SetOrigin(glm::vec3(32.5, 125.5, 43.7));
    //fluidsph.SetOrigin(glm::vec3(32.5, 160.5, 43.7)); //lena
    
    //camera.PlaceTo(glm::vec3(37.366, 128.401, 41.44)); //video 1
    camera.PlaceTo(glm::vec3(29.798, 126.626, 43.954)); //video 2
    //camera.PlaceTo(glm::vec3(35.673, 124.497, 39.898)); //video 3
    fluidsph.Initialize(1000);


    glEnable(GL_DEPTH_TEST);
    glReadBuffer(GL_BACK);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    UIinit(window);

    unsigned char* buff = new unsigned char[SCREEN_HEIGHT * SCREEN_WIDTH * 3];
    char filename[50];
    unsigned int framenum = 0;
    std::stringstream ss_filename;

    //render loop
    while (!glfwWindowShouldClose(window))
    {
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        processInput(window);
        glClearColor(0.3, 0.3, 0.3, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        UIbegin();
        
        if (!disabled)
        {
            ImGui::Begin("New window");

            ImGui::InputFloat3("Camera pos:", (float*)&camera.Position);

            ImGui::ColorEdit3("clear color", (float*)&clear_color);
            grid.SetColor(clear_color);
            ImGui::InputFloat3("dimensions", (float*)&dimensions);
            if (ImGui::Button("Update grid"))
            {
                grid.UpdateGrid(dimensions);
            }
            ImGui::InputText("picture path", picture_path, 100);
            if (ImGui::Button("Load picture"))
            {
                grid.UpdateHeightMap(picture_path);
                grid.UpdateGrid(dimensions);
            }
            ImGui::Checkbox("render boundary", fluidsph.GetRenderBoundary());

            ImGui::NewLine();

            ImGui::InputFloat("Friction", fluidsph.GetMu(), 0.01f);
            ImGui::InputFloat("Ks", fluidsph.GetKs(), 1000.0f);
            ImGui::InputFloat("Kd", fluidsph.GetKd(), 1000.0f);
            ImGui::InputFloat("Mass", fluidsph.GetMass(), 0.001f);
            ImGui::InputFloat("Visc", fluidsph.GetVisc(), 0.001f);
            ImGui::InputFloat("SurfTens", fluidsph.GetSurfTen(), 0.0001f, 0.0f, "%.4f");
            ImGui::InputFloat("p0", fluidsph.Getp0(), 1.0f);
            ImGui::InputFloat("damping", fluidsph.GetDamping(), 0.1f);
            float g[3];
            glm::vec3* gr = fluidsph.GetGrav();
            g[0] = gr->x;
            g[1] = gr->y;
            g[2] = gr->z;
            ImGui::InputFloat3("Gravity", g);
            gr->x = g[0];
            gr->y = g[1];
            gr->z = g[2];
            ImGui::NewLine();
            ImGui::InputInt("particle_id", &g_part_id);
            FluidParticle fp = fluidsph.GetParticle(g_part_id);
            ImGui::InputFloat3("fBoundary", (float*)&fp.fBoundary);
            ImGui::InputFloat("d", &fp.shortest);
            ImGui::InputFloat3("Pos", (float*)&fp.Position);
            ImGui::InputFloat3("Vel", (float*)&fp.Velocity);
            ImGui::InputFloat3("Acc", (float*)&fp.Acceleration);
            ImGui::InputFloat("Dens", &fp.Density);
            ImGui::InputFloat("Pres", &fp.Pressure);
            ImGui::InputFloat3("PresF", (float*)&fp.PressureForce);
            //ImGui::InputFloat3("fPressTmp", (float*)&fp.fPressTmp);
            ImGui::InputInt("NeighbId", &fp.NeighbId, 0);
            ImGui::InputFloat3("GravF", (float*)&fp.GravityForce);
            ImGui::InputFloat3("SurfF", (float*)&fp.SurfaceForce);
            ImGui::InputFloat3("SurfNorm", (float*)&fp.SurfaceNormal);
            ImGui::InputFloat3("ViscF", (float*)&fp.ViscosityForce);


            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
            ImGui::End();
        }

        
        
        fluidsph.Run(grid);// simulation

        glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)SCREEN_WIDTH / SCREEN_HEIGHT, 0.01f, 1000.0f);
        glm::mat4 view = camera.GetViewMatrix();
        shader.use();
        shader.SetProjection(projection);
        shader.SetView(view);
        shader.setVec3("viewerPos", camera.Position);

        grid.Draw(shader);
        fluidsph.Draw(shader, g_part_id);
        
        UIend();

        glReadPixels(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT, GL_BGR, GL_UNSIGNED_BYTE, buff);
        ss_filename << "render/frame_" << std::to_string(framenum) << ".bmp";
        framenum++;
        //pixelsToBmp(ss_filename.str().c_str(), buff);
        ss_filename.str("");

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    
    UIshutdown();
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

    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS && disabled)
    {
        fluidsph.AddParticles(125);
        std::cout << "Added 125 particles\n";
    }

    if (button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_PRESS)
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

GLFWwindow* init()
{
    if (!glfwInit())
        exit(-1);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Erosion", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create glfw window\n";
        glfwTerminate();
        exit(-1);
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
        exit(-1);
    }

    return window;
}

void UIinit(GLFWwindow* window)
{
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
}

void UIshutdown()
{
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

void UIbegin()
{
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

void UIend()
{
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void pixelsToBmp(const char* filename, const unsigned char* pixels)
{
#pragma warning(suppress : 4996)
    FILE* out = fopen(filename, "wb");
    if (!out)
    {
        fprintf(stderr, "can't open file %s\n", filename);
        exit(1);
    }
    BITMAPFILEHEADER bitmapFileHeader;
    BITMAPINFOHEADER bitmapInfoHeader;

    bitmapFileHeader.bfType = 0x4D42;
    bitmapFileHeader.bfSize = SCREEN_WIDTH * SCREEN_HEIGHT * 3;
    bitmapFileHeader.bfReserved1 = 0;
    bitmapFileHeader.bfReserved2 = 0;
    bitmapFileHeader.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);

    bitmapInfoHeader.biSize = sizeof(BITMAPINFOHEADER);
    bitmapInfoHeader.biWidth = SCREEN_WIDTH - 1;
    bitmapInfoHeader.biHeight = SCREEN_HEIGHT - 1;
    bitmapInfoHeader.biPlanes = 1;
    bitmapInfoHeader.biBitCount = 24;
    bitmapInfoHeader.biCompression = BI_RGB;
    bitmapInfoHeader.biSizeImage = 0;
    bitmapInfoHeader.biXPelsPerMeter = 0; // ?
    bitmapInfoHeader.biYPelsPerMeter = 0; // ?
    bitmapInfoHeader.biClrUsed = 0;
    bitmapInfoHeader.biClrImportant = 0;

    fwrite(&bitmapFileHeader, sizeof(BITMAPFILEHEADER), 1, out);
    fwrite(&bitmapInfoHeader, sizeof(BITMAPINFOHEADER), 1, out);
    fwrite(pixels, SCREEN_WIDTH * SCREEN_HEIGHT * 3, 1, out);
    fclose(out);
}
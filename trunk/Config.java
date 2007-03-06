import java.awt.event.KeyEvent;

/**
 * A <code>Config</code> object defines a number of constants used
 * by the game.
 * Editing this class could lead to serious gaming problems.
 */
public class Config
{
    /*
     * Key constants
     */

    // Move right key
    public static final int MOVE_RIGHT = KeyEvent.VK_RIGHT;

    // Move left key
    public static final int MOVE_LEFT = KeyEvent.VK_LEFT;

    // Move down key
    public static final int MOVE_DOWN = KeyEvent.VK_DOWN;

    // Move up key
    public static final int MOVE_UP = KeyEvent.VK_UP;

    /*
     * Game constants
     */

    // The width of the tiles
    public static final int GRID_WIDTH = 30;

    // The height of the tiles
    public static final int GRID_HEIGHT = 30;

    // The width of a game
    public static final int GAME_WIDTH = 15 * GRID_WIDTH;

    // The height of a game
    public static final int GAME_HEIGHT = 11 * GRID_HEIGHT;

    // The maximum number of players
    public static final int MAX_PLAYERS = 8;


    /*
     * Server constants
     */

    // The unique id given to the server. Clients can
    // not have this id!
    public static final int SERVER_ID = 0;

    // The port to start a server on
    public static final int SERVER_PORT = 4444;

    // Path to save log files to
    public static final String LOG_PATH = "LOGS/";

    // Server hostname
    public static final String SERVER_ADDR = "localhost";


    /*
     * Client constants
     */

    // The pad where scheme files can be found
    public static final String SCHEME_PATH = "DATA/SCHEMES/";


    /*
     * Command constants
     */

    // The id of an InitCommand
    public static final int CMD_INIT = 1;

    // The id of a NewPlayerCommand
    public static final int CMD_NEW = 2;

    // The id of a GameCommand
    public static final int CMD_GAME = 3;

    // The id of a PlayerCommand
    public static final int CMD_PLAYER = 4;

    /**
     * Action constants
     */

    // Action move right
    public static final int ACT_MOVE_RIGHT = 1;

    // Action move left
    public static final int ACT_MOVE_LEFT = 2;

    // Action move down
    public static final int ACT_MOVE_DOWN = 3;

    // Action move up
    public static final int ACT_MOVE_UP = 4;

    // Action nothing
    public static final int ACT_NOTHING = 5;

    // Action leave
    public static final int ACT_LEAVE = 6;
}
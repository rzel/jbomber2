import java.io.*;

/**
 * A <code>NewPlayerCommand</code> is a command that is send by a client
 * whenever an <code>InitCommand has been send by the server, or whenever a
 * client has recieved another <code>NewPlayerCommand</code>.
 * This command contains player specific data, which means the position of a
 * player, the nickname of a player and a flag that indicates whether the
 * player had previously joinded the game.
 */
public class NewPlayerCommand extends Command
{
    // Something
    public static final long serialVersionUID = 1;


    // The horizontal position of the player
    private int x_pos;

    // The vertical position of the player
    private int y_pos;

    // The nickname of the player
    private String nickname;

    // true indicates a new player
    private boolean new_player;

    /**
     * Creates a <code>NewPlayerCommand</code> with the position of
     * a player, the nickname of a player and a <code>boolean</code> flag.
     * When <code>true</code> when this flag indicates that the client
     * that send this command was not previously connected (the command
     * was send DIRECTLY after an <code>InitCommand</code>).
     * When <code>False</code> this flag indicates that the client had
     * already joined the game.
     *
     * @param sender_id  The unique id of the player who send this command.
     * @param x_pos      The horizontal position of the player.
     * @param y_pos      The vertical position of the player.
     * @param nickname   The nickname of the player who send this command.
     * @param new_player <code>true</code> indicates a new player.
     */
    public NewPlayerCommand(int sender_id, int x_pos, int y_pos,
                            String nickname, boolean new_player)
    {
        super(Config.CMD_NEW, sender_id);
        this.x_pos      = x_pos;
        this.y_pos      = y_pos;
        this.nickname   = nickname;
        this.new_player = new_player;
    }

    /**
     * Returns the horizontal position of the player.
     *
     * @return The horizontal position if the player.
     */
    public int getXPos()
    {
        return x_pos;
    }

    /**
     * Returns the vertical position of the player.
     *
     * @return The vertical position if the player.
     */
    public int getYPos()
    {
        return y_pos;
    }

    /**
     * Returns the nickname of the player.
     *
     * @return The nickname if the player.
     */
    public String getNickname()
    {
        return nickname;
    }

    /**
     * Returns whether this was a new player.
     *
     * @return <code>true</code> when this player was new,
     *         otherwise <code>false</code>.
     */
    public boolean getNewPlayer()
    {
        return new_player;
    }
}